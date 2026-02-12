library(ArchR)
library(Seurat)
library(visutils)
source('actions/ovary/bin/archr_utils.R')

addArchRGenome("hg38")
addArchRThreads(threads = 1)
# project = ArchR::loadArchRProject('work/archr2/02_annotate_out')
gsm = readRDS('work/archr2/02_annotate_out/GeneScoreMatrix.rds')
samples = read.csv('actions/samples.csv')
# samples$arrows = normalizePath(paste0('work/archr/01_arrows/arrows/',samples$sample_id,'/',samples$sample_id,'.arrow'))
# file.exists(samples$arrows)
# gene.descr = readRDS('/nfs/cellgeni/pasham/code/gene.descr.rds')
# 
# outdir = 'work/archr2/02_annotate'

# annotate cells based on Valentinas approach
# https://github.com/ventolab/Human-ReproductiveTract-Development-Atlas/blob/main/preprocessing/scatacseq/SingleSample_Analysis_ATAC_ArchR.ipynb

# geneAnnotation = myCreateGeneAnnotation('hg38','/software/cellgen/cellgeni/refdata_10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf')
# saveRDS(geneAnnotation,'geneAnnotation_2020A.rds')
geneAnnotation = readRDS('geneAnnotation_2020A.rds')

project <- ArchRProject(
  ArrowFiles = samples$arrows,
  outputDirectory = outdir,
  copyArrows = TRUE
)

project

# with links matrix it too large
geneAnnotation$genes = geneAnnotation$genes[geneAnnotation$genes$gene_type == 'protein_coding' | geneAnnotation$genes$symbol %in%  ArchR::getGenes()$symbol, ]
length(geneAnnotation$genes) # 36559 -> 22073

project = addGeneScoreMatrix(project,genes = geneAnnotation$genes,force=TRUE)
#saveArchRProject(project,outputDirectory = 'work/archr2/02_annotate_out')
# add metadata ######
samples_info = read.csv('data/Table1-dataset_metadata - scATACseq_SangerPediatric.csv',row.names = 1)
project$Age = samples_info[project$Sample,'Age']
project$Donor = samples_info[project$Sample,'Donor']


# QC plots #############
df <- getCellColData(project, select = c("log10(nFrags)", "TSSEnrichment"))
p_raw <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 10, lty = "dashed") + geom_vline(xintercept = log10(3500), lty = "dashed")
ggsave('actions/ovary/figures/arch_qc/nfrags2tss.pdf',p_raw)


p_tss <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)

ggsave('actions/ovary/figures/arch_qc/tsse.pdf',p_tss)


p_tss_v <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
ggsave('actions/ovary/figures/arch_qc/tsse_violin.pdf',p_tss_v)




p_frags <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)

ggsave('actions/ovary/figures/arch_qc/nfrags.pdf',p_frags)


p_frag_2 <- plotFragmentSizes(ArchRProj = project)
ggsave('actions/ovary/figures/arch_qc/frags_sizes.pdf',p_frag_2)


p_tss_2 <- plotTSSEnrichment(ArchRProj = project)
ggsave('actions/ovary/figures/arch_qc/tss_enich.pdf',p_tss_2)


# dim reductions/clustering/umap #############
project <- addIterativeLSI(
  ArchRProj = project,
  useMatrix = "TileMatrix", 
  name = "TileMatrix_LSI", 
  iterations = 2, 
  clusterParams = list(resolution = c(2), sampleCells = 10000, maxClusters = 6, n.start = 10), 
  varFeatures = 25000, 
  dimsToUse = 1:30, 
  LSIMethod = 2
)


project = addHarmony(
  ArchRProj = project,
  reducedDims = "TileMatrix_LSI",
  name = "TileMatrix_LSI_harmony",
  groupBy = "Donor"
)

for(dim in 'TileMatrix_LSI'){#names(project@reducedDims)){
  print(dim)
  for(r in c(4,6,8,10)){#,1,2)){
    print(r)
    project <- addClusters(
      input = project,
      reducedDims = dim,
      method = "Seurat",
      name = paste0(dim,'_cl_',r),
      resolution = r, 
      #maxClusters = 25, 
      knnAssign = 20, 
      force = TRUE
    )
  }
  
  # project <- addUMAP(
  #   ArchRProj = project, 
  #   reducedDims = dim, 
  #   name = paste0(dim,'_UMAP'), 
  #   nNeighbors = 30, 
  #   minDist = 0.5, 
  #   metric = "cosine"
  # )
}

# integrate ############
rna = schard::h5ad2seurat('data/ref_rawcounts.h5ad')
y=rna@assays$RNA@counts['LHX2',] 
boxplot(log1p(y/rna$nCount_RNA*1e4) ~ splitSub(rna$coarse_annotation,'_',1),las=2)
sort(tapply(log1p(y/rna$nCount_RNA*1e4) , splitSub(rna$coarse_annotation,'_',1),mean))
sort(tapply(y>0 , splitSub(rna$coarse_annotation,'_',1),mean))


project <- addGeneIntegrationMatrix(
  ArchRProj = project, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "TileMatrix_LSI_harmony",
  seRNA = rna,
  addToArrow = FALSE,
  force = TRUE,
  groupRNA = 'coarse_annotation',
  nameCell = "coarse_annotation_archr_h_cell",
  nameGroup = "coarse_annotation_archr_h",
  nameScore = "coarse_annotation_archr_h_score"
)


project <- addGeneIntegrationMatrix(
  ArchRProj = project, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "TileMatrix_LSI",
  seRNA = rna,
  addToArrow = FALSE,
  force = TRUE,
  groupRNA = 'coarse_annotation',
  nameCell = "coarse_annotation_archr_cell",
  nameGroup = "coarse_annotation_archr",
  nameScore = "coarse_annotation_archr_score"
)

obs = schard::h5ad2data.frame('data/ref_rawcounts.h5ad','obs')
# these are not really harmonized
celltypes = unique(obs[,c('lineage','coarse_annotation')])
celltypes$class = visutils::splitSub(celltypes$coarse_annotation,'_',1)
table(celltypes$class,celltypes$lineage)
celltypes[celltypes$class=='PV',]

project$class   = visutils::splitSub(project$coarse_annotation_archr,'_',1)
project$class_h = visutils::splitSub(project$coarse_annotation_archr_h,'_',1)
project$class[project$class %in% c('Theca','Stroma','FibC3')] = 'Mesenchymal'
project$class_h[project$class_h %in% c('Theca','Stroma','FibC3')] = 'Mesenchymal'

# celltypes missed in ATAC
setdiff(celltypes$class,visutils::splitSub(project$coarse_annotation_archr,'_',1))
setdiff(celltypes$coarse_annotation,project$coarse_annotation_archr)

# umaps ##########################
t = table(project$class,project$TileMatrix_LSI_cl_1)
cl2class = setNames(rownames(t)[apply(t,2,which.max)],colnames(t))
project$class_clean = cl2class[project$TileMatrix_LSI_cl_1]


png('work/archr2/02_annotate_out/figures/01_whole_TileMatrix_LSI_UMAP.png',width = 4.5*3,h=12,units = 'in',res = 300)
par(mfrow=c(3,3),mar=c(1,1,1,6),oma=c(0,0,1,0))
umap = project@embeddings$TileMatrix_LSI_UMAP$df
visutils::plotVisium(umap,project$Donor,t='xy',cex=0.2,plot.legend = T)
visutils::plotVisium(umap,project$DoubletScore,t='xy',cex=0.2,plot.legend = T,order.points.by.z = T)
visutils::plotVisium(umap,project$class,t='xy',cex=0.1,label.clusters = T,randomize.points = T)
visutils::plotVisium(umap,project$coarse_annotation_archr,t='xy',cex=0.1,label.clusters = T,randomize.points = T,legend.args = list(ncol=2))
plot.new()
visutils::plotVisium(umap,project$TileMatrix_LSI_cl_1,t='xy',cex=0.1,label.clusters = T,randomize.points = T)
visutils::plotVisium(umap,project$class_clean,t='xy',cex=0.1,label.clusters = T,randomize.points = T,main='Cleaned by TileMatrix_LSI_cl_1 classes')
mtext("TileMatrix_LSI_UMAP",side = 3,outer = TRUE)
dev.off()

png('work/archr2/02_annotate_out/figures/01_whole_TileMatrix_LSI_harmony_UMAP.png',width = 4.5*3,h=8,units = 'in',res = 300)
par(mfrow=c(2,3),mar=c(1,1,1,6),oma=c(0,0,1,0))
umap = project@embeddings$TileMatrix_LSI_harmony_UMAP$df
visutils::plotVisium(umap,project$Donor,t='xy',cex=0.2,plot.legend = T)
visutils::plotVisium(umap,project$DoubletScore,t='xy',cex=0.2,plot.legend = T,order.points.by.z = T)
visutils::plotVisium(umap,project$class_h,t='xy',cex=0.1,label.clusters = T,randomize.points = T)
visutils::plotVisium(umap,project$coarse_annotation_archr_h,t='xy',cex=0.1,label.clusters = T,randomize.points = T,legend.args = list(ncol=2))
plot.new()
visutils::plotVisium(umap,project$TileMatrix_LSI_harmony_cl_1,t='xy',cex=0.1,label.clusters = T,randomize.points = T)
mtext("TileMatrix_LSI_harmony_UMAP",side = 3,outer = TRUE)
dev.off()

saveArchRProject(project,outputDirectory = 'work/archr2/02_annotate_out')

# harmony splits granulosa into parts... so I'll not use it
# rename clusters bu classes
umap = project@embeddings$TileMatrix_LSI_UMAP$df
visutils::plotVisium(umap,project$class_clean,t='xy',cex=0.4,label.clusters = T,randomize.points = T,show.cluster.sizes = T)
visutils::plotVisium(umap,project$Donor,t='xy',cex=0.4,label.clusters = F,randomize.points = T,legend.args = list(title='Donor'))
visutils::plotVisium(umap,project$class_h=='PV',t='xy',order.points.by.z = T,cex=(project$class_h=='PV')+0.1)
visutils::plotVisium(umap,project$TileMatrix_LSI_cl_1=='C3',t='xy',order.points.by.z = T,cex=(project$TileMatrix_LSI_cl_1=='C3')+0.1)


markerGenes = c('LHX2','UPK3B','MSLN', 'LRRN4', 'KLK11')

gs = ArchR::getMatrixFromProject(project,useMatrix = 'GeneScoreMatrix')
saveRDS(gs,'work/archr2/02_annotate_out/GeneScoreMatrix.rds')
mtx = assay(gs,'GeneScoreMatrix')
rownames(mtx) = rowData(gs)$name
colnames(mtx) = colnames(gs)
mtx = mtx[,project$cellNames]
y = mtx['KLK11',]
boxplot(log10(1+y) ~ project$class_clean)

plotVisium(umap,y,t='xy',order.points.by.z = T,label.clusters = project$class_clean,zfun = log1p,cex=0.2)

markerGenes %in% rownames(mtx)

p <- plotEmbedding(
  ArchRProj = project, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "TileMatrix_LSI_UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
p



# annotated per coarse celltypes groups #########
project = ArchR::loadArchRProject('work/archr2/02_annotate_out')
prj = ArchR::loadArchRProject('work/archr_with_fetal/02_annotate/')
table(prj$Sample[prj$potential_OSE])
table(prj$Donor[prj$potential_OSE])

table(project$TileMatrix_LSI_cl_10,project$cellNames %in% prj$cellNames[prj$potential_OSE]) # C13

table(project$class_clean[project$TileMatrix_LSI_cl_10=='C13']) # Granulosa

project$class_clean[project$TileMatrix_LSI_cl_10=='C13'] = 'OSE'
table(project$class_clean)

path = 'work/archr2/03_per_class_annotation'
rna = schard::h5ad2seurat('data/ref_rawcounts.h5ad')
rna$class = visutils::splitSub(rna$coarse_annotation,'_',1)
rna$class[rna$class %in% c('Theca','Stroma','FibC3')] = 'Mesenchymal'

#dir.create(path)
setdiff(project$class_clean,rna$class)
ArchR::addArchRVerbose(FALSE)

for(class in "Granulosa"){# unique(project$class_clean)){
  print(class)
  gc()
  subset = subsetArchRProject(
    ArchRProj = project,
    cells = project$cellNames[project$class_clean==class],
    outputDirectory = paste0(path,"/",class),
    dropCells = TRUE,
    force = TRUE
  )
  
  subset = addIterativeLSI(
    ArchRProj = subset,
    useMatrix = "TileMatrix", 
    name = "TileMatrix_LSI", 
    iterations = 2, 
    clusterParams = list(resolution = c(2), sampleCells = 10000, maxClusters = 6, n.start = 10), 
    varFeatures = 25000, 
    dimsToUse = 1:30, 
    LSIMethod = 2,
    force = TRUE
  )
  
  
  subset = addHarmony(
    ArchRProj = subset,
    reducedDims = "TileMatrix_LSI",
    name = "TileMatrix_LSI_harmony",
    groupBy = "Donor",
    force = TRUE
  )
  
  for(dim in names(subset@reducedDims)){
    print(dim)
    for(r in c(1,2)){
      print(r)
      subset <- addClusters(
        input = subset,
        reducedDims = dim,
        method = "Seurat",
        name = paste0(dim,'_cl_',r),
        resolution = r, 
        #maxClusters = 25, 
        knnAssign = 20, 
        force = TRUE
      )
    }
    
    subset <- addUMAP(
      ArchRProj = subset,
      reducedDims = dim,
      name = paste0(dim,'_UMAP'),
      nNeighbors = 30,
      minDist = 0.5,
      metric = "cosine",
      force = TRUE
    )
  }
  
  rna_sub = rna[,rna$class==class]
  
  subset <- addGeneIntegrationMatrix(
    ArchRProj = subset, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "TileMatrix_LSI_harmony",
    seRNA = rna_sub,
    addToArrow = FALSE,
    force = TRUE,
    groupRNA = 'coarse_annotation',
    nameCell = "class-coarse_annotation_archr_h_cell",
    nameGroup = "class-coarse_annotation_archr_h",
    nameScore = "class-coarse_annotation_archr_h_score"
  )
  
  
  subset <- addGeneIntegrationMatrix(
    ArchRProj = subset, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "TileMatrix_LSI",
    seRNA = rna_sub,
    addToArrow = FALSE,
    force = TRUE,
    groupRNA = 'coarse_annotation',
    nameCell = "class-coarse_annotation_archr_cell",
    nameGroup = "class-coarse_annotation_archr",
    nameScore = "class-coarse_annotation_archr_score"
  )
  
  
  saveArchRProject(subset,outputDirectory = paste0(path,"/",class))
}

# check classes #############
gs = readRDS('work/archr_with_fetal_atac/02_annotate/GeneScoreMatrix.rds')

# _add to major object #######
project$coarse_annotation_archr_clean = NA
project$coarse_annotation_archr_h_clean = NA

for(cl in unique(project$class_clean)){
  print(cl)
  if(cl == 'OSE'){
    f = project$class_clean == 'OSE'
    project@cellColData[f,'coarse_annotation_archr_h_clean'] = 'OSE'
    project@cellColData[f,'coarse_annotation_archr_clean'] = 'OSE'
    next
  }
  clp = ArchR::loadArchRProject(paste0('work/archr2/03_per_class_annotation/',cl))
  clp$class.coarse_annotation_archr_clean = getMajorCelltype(clp$TileMatrix_LSI_cl_2,clp$class.coarse_annotation_archr)
  clp$class.coarse_annotation_archr_h_clean = getMajorCelltype(clp$TileMatrix_LSI_harmony_cl_2,clp$class.coarse_annotation_archr_h)
  project@cellColData[clp$cellNames,'coarse_annotation_archr_h_clean'] = clp$class.coarse_annotation_archr_h_clean
  project@cellColData[clp$cellNames,'coarse_annotation_archr_clean'] = clp$class.coarse_annotation_archr_clean
}

(t=table(project$coarse_annotation_archr_clean,useNA='always'))
(t=table(project$coarse_annotation_archr_h,useNA='always'))


plotVisium(project@embeddings$TileMatrix_LSI_UMAP$df,project$coarse_annotation_archr_clean,cex=0.3,label.clusters = T,t='xy')
plotVisium(project@embeddings$TileMatrix_LSI_harmony_UMAP$df,project$coarse_annotation_archr_h_clean,cex=0.3,label.clusters = T,t='xy')

saveArchRProject(project,outputDirectory = 'work/archr2/02_annotate_out')

# _granulosa ############
gra = ArchR::loadArchRProject('work/archr2/03_per_class_annotation/Granulosa/')


grags = assay(gs,'GeneScoreMatrix')[,gra$cellNames]
rownames(grags) = rowData(gs)$name
colnames(grags) = gra$cellNames

lhx2 = grags['LHX2',]
amh = grags['AMH',]

par(mar=c(1,1,1,15),xpd=NA)
plotVisium(gra@embeddings$TileMatrix_LSI_harmony_UMAP$df,gra$class.coarse_annotation_archr_h,t='xy',label.clusters = T,show.cluster.sizes = T,cex=1)
plotVisium(gra@embeddings$TileMatrix_LSI_harmony_UMAP$df,gra$coarse_annotation_archr_h,t='xy',label.clusters = T,show.cluster.sizes = T,cex=1)
plotVisium(gra@embeddings$TileMatrix_LSI_harmony_UMAP$df,gra$DoubletScore,t='xy',label.clusters = T,show.cluster.sizes = T,cex=0.5)
plotVisium(gra@embeddings$TileMatrix_LSI_harmony_UMAP$df,gra$Donor,t='xy',label.clusters = T,show.cluster.sizes = T,cex=0.5)
plotVisium(gra@embeddings$TileMatrix_LSI_harmony_UMAP$df,gra$TileMatrix_LSI_harmony_cl_2,t='xy',label.clusters = T,show.cluster.sizes = T,cex=0.5)

plotVisium(gra@embeddings$TileMatrix_LSI_UMAP$df,gra$class.coarse_annotation_archr,t='xy',label.clusters = T,show.cluster.sizes = T,cex=1)
plotVisium(gra@embeddings$TileMatrix_LSI_UMAP$df,gra$coarse_annotation_archr,t='xy',label.clusters = T,show.cluster.sizes = T,cex=1)
plotVisium(gra@embeddings$TileMatrix_LSI_UMAP$df,gra$DoubletScore,t='xy',label.clusters = T,show.cluster.sizes = T,cex=0.5)
plotVisium(gra@embeddings$TileMatrix_LSI_UMAP$df,gra$Donor,t='xy',label.clusters = T,show.cluster.sizes = T,cex=0.5)

# _endothelial ############    
# here we saw some potential OSE

endo = ArchR::loadArchRProject('work/archr2/03_per_class_annotation/Endo')

endogs = assay(gs,'GeneScoreMatrix')[,endo$cellNames]
rownames(endogs) = rowData(gs)$name
colnames(endogs) = endo$cellNames

lhx2 = endogs['LHX2',]


par(mar=c(1,1,1,15),xpd=NA)
plotVisium(endo@embeddings$TileMatrix_LSI_harmony_UMAP$df,endo$class.coarse_annotation_archr_h,t='xy',label.clusters = T,show.cluster.sizes = T,cex=0.5)
plotVisium(endo@embeddings$TileMatrix_LSI_harmony_UMAP$df,endo$coarse_annotation_archr_h,t='xy',label.clusters = T,show.cluster.sizes = T,cex=0.5)
plotVisium(endo@embeddings$TileMatrix_LSI_harmony_UMAP$df,endo$Donor,t='xy',label.clusters = T,show.cluster.sizes = T,cex=0.5)
plotVisium(endo@embeddings$TileMatrix_LSI_harmony_UMAP$df,lhx2,t='xy',cex=0.5)

plotVisium(endo@embeddings$TileMatrix_LSI_UMAP$df,endo$class.coarse_annotation_archr,t='xy',label.clusters = T,show.cluster.sizes = T,cex=0.5)
plotVisium(endo@embeddings$TileMatrix_LSI_UMAP$df,endo$coarse_annotation_archr,t='xy',label.clusters = T,show.cluster.sizes = T,cex=0.5)
plotVisium(endo@embeddings$TileMatrix_LSI_UMAP$df,endo$Donor,t='xy',label.clusters = T,show.cluster.sizes = T,cex=0.5)
plotVisium(endo@embeddings$TileMatrix_LSI_UMAP$df,lhx2,t='xy',label.clusters = T,show.cluster.sizes = T,cex=0.5)

table(endo$class.coarse_annotation_archr_h,endo$class.coarse_annotation_archr)

# _Immune ############
# here we saw some potential OSE

imm = ArchR::loadArchRProject('work/archr2/03_per_class_annotation/Immune/')

immgs = assay(gs,'GeneScoreMatrix')[,imm$cellNames]
rownames(immgs) = rowData(gs)$name
colnames(immgs) = imm$cellNames

lhx2 = immgs['LHX2',]


par(mar=c(1,1,1,15),xpd=NA)
plotVisium(imm@embeddings$TileMatrix_LSI_harmony_UMAP$df,imm$class.coarse_annotation_archr_h,t='xy',label.clusters = T,show.cluster.sizes = T,cex=1)
plotVisium(imm@embeddings$TileMatrix_LSI_harmony_UMAP$df,imm$coarse_annotation_archr_h,t='xy',label.clusters = T,show.cluster.sizes = T,cex=1)
plotVisium(imm@embeddings$TileMatrix_LSI_harmony_UMAP$df,imm$Donor,t='xy',label.clusters = T,show.cluster.sizes = T,cex=0.5)

plotVisium(imm@embeddings$TileMatrix_LSI_UMAP$df,imm$class.coarse_annotation_archr,t='xy',label.clusters = T,show.cluster.sizes = T,cex=0.5)
plotVisium(imm@embeddings$TileMatrix_LSI_UMAP$df,imm$coarse_annotation_archr,t='xy',label.clusters = T,show.cluster.sizes = T,cex=0.5)
plotVisium(imm@embeddings$TileMatrix_LSI_UMAP$df,imm$Donor,t='xy',label.clusters = T,show.cluster.sizes = T,cex=0.5)

imageWithText(table(imm$class.coarse_annotation_archr_h,imm$class.coarse_annotation_archr))

# _Oocytes ##############
rna_obs = schard::h5ad2data.frame('data/ref_rawcounts.h5ad','/obs')
rnact = table(rna_obs$coarse_annotation)
oocytesn = rnact[grep('oo',names(rnact),ignore.case = T)]
oocytes = names(oocytesn)
table(project$coarse_annotation_archr)[oocytes]

oo = (project$class=='Germcell') + 0
o = order(oo,decreasing = F)
par(mar=c(0,0,1,14))
plotVisium(project@embeddings$TileMatrix_LSI_UMAP$df[o,],project$class[o],t='xy',cex=oo[o]+0.2,label.clusters = T)
#plotVisium(project@embeddings$TileMatrix_LSI_UMAP$df,oo,t='xy',cex=oo+0.2,order.points.by.z = T)


oo = (project$class_h=='Germcell') + 0
o = order(oo,decreasing = F)
par(mar=c(0,0,1,14))
plotVisium(project@embeddings$PeakMatrix_LSI_harmony_UMAP$df[o,],project$class_h[o],t='xy',cex=oo[o]+0.2,label.clusters = T)

assayNames(gsm)
m = unwrapAssay(gsm)
f = gsm$class_clean
f[gsm$class =='Germcell'] = 'Germcell'
f = as.factor(f)

gids=c('DAZL', 'DDX4', 'FIGLA', 'GDF9', 'ZP3', 'NOBOX')
setdiff(gids,rownames(m))
par(mar=c(8,6,1,10),bty='n',las=2)
plotExpAsDots(m,f,gids,max.cex=4)
plotExpAsDots(m,gsm$TileMatrix_LSI_cl_10,gids,max.cex=4,main='TileMatrix_LSI_cl_10')
table(gsm$TileMatrix_LSI_cl_10,f)['C2',]

project = addModuleScore(project,useMatrix = 'GeneScoreMatrix',name='module',features = list(germcell=gids))

plotVisium(project@embeddings$TileMatrix_LSI_UMAP$df,project$module.germcell,t='xy',cex=0.3,
           label.clusters = project$class,order.points.by.z = T)

z = project$module.germcell>=sort(project$module.germcell,decreasing = T)[1000]
t = table(project$class,z)
sort(t[,2]/rowSums(t))

plotVisium(project@embeddings$TileMatrix_LSI_UMAP$df,z,t='xy',cex=0.3,
           label.clusters = project$class,order.points.by.z = T)
par(mar=c(0,0,1,6),bty='n',las=2)
plotVisium(project@embeddings$TileMatrix_LSI_UMAP$df,project$module.germcell,t='xy',cex=0.3,
           label.clusters = project$class,order.points.by.z = T,zfun=log1p,main='module.germcell')
boxplot(project$module.germcell ~ project$class)
boxplot(project$module.germcell ~ project$TileMatrix_LSI_cl_10)

# marker genes for classes #######
# markersGSclass <- getMarkerFeatures(
#   ArchRProj = project,
#   useMatrix = "GeneScoreMatrix",
#   groupBy = "class_clean",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# saveRDS(markersGSclass,'work/archr2/02_annotate_out/markersGSclass.rds')

markersGSclass = readRDS('work/archr2/02_annotate_out/markersGSclass.rds')
rownames(markersGSclass) = make.unique(rowData(markersGSclass)$name)


top5 = myGetMarkers(markersGSclass, cutOff = "FDR <= 0.05 & Log2FC >= 0.6",n=6)
table(top5$comparison)

top5l = split(top5$name,top5$comparison)
top5l$OSE = unique(c(top5l$OSE,'UPK3B','LHX2','LRRN4'))

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGSclass[unlist(top5l),], 
  cutOff = "FDR <= 10", 
  transpose = TRUE,clusterCols=FALSE,nLabel = 1000
)
heatmapGS

gids = unlist(top5l)
genect = unlist(lapply(names(top5l),function(n)rep(n,length(top5l[[n]]))))
setdiff(gids,rownames(m))

ctcols = char2col(f)
plotExpAsDots(m,f,gids=unlist(top5l),scale_exp_mar=1,max.cex = 2,rowAnnWidth=0.4,
              colColours=ctcols[levels(f)],
              rowColours=ctcols[genect])



umap=project@embeddings$PeakMatrix_LSI_UMAP$df
umap=project@embeddings$TileMatrix_LSI_UMAP$df
par(mar=c(0,0,1,22))
plotVisium(umap,project$class_clean,t='xy',cex=0.2,label.clusters = T,show.cluster.sizes = T)
plotVisium(umap,project$coarse_annotation_archr_clean,t='xy',cex=0.4,label.clusters = T,show.cluster.sizes = T,randomize.points = T)

# _gene from Luz ###########
class_markers = list('Germ cells' = c('DAZL', 'ZP3'), 
     'Granulosa' = c("AMH", 'RDH10', 'FOXL2', 'KITLG'),
     'Mesothelial' = c( 'LHX2', 'LRRN4', 'UPK3B'),  
     'Mesenchymal' = c('DCN', 'COL1A2', 'PDGFRA'),
     'Theca' = c("LHCGR", "PTCH1", 'INSL3', 'THBD', 'WDFY1'), 
     'PV' = c('RGS5',  'RERGL', 'KCNAB1', 'MYH11'),
     'Endothelial' = c('PECAM1', 'CDH5'))
class_markerso = names(class_markers)
class_markers = unlist(class_markers)
class_markers = data.frame(class=names(class_markers),marker=class_markers)
class_markers$class = sub("\\d$",'',class_markers$class)


setdiff(unlist(class_markers),rownames(markersGSclass))
cols=list(lineage=char2col(class_markerso)[class_markerso])

linord = c("Granulosa","OSE" ,"Mesenchymal","PV","Endo","Immune","Neural")
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGSclass[class_markers$marker,linord], 
  cutOff = "FDR <= 10", 
  transpose = TRUE,clusterCols=FALSE,nLabel = 1000
)
heatmapGS = HeatmapAnnotation(lineage=class_markers$class,col=cols,which='column')  %v% heatmapGS  

heatmapGS 

project$coarse_annotation_archr[grep('oo',project$coarse_annotation_archr,ignore.case = T)]
project$class[grep('oo',project$coarse_annotation_archr,ignore.case = T)]
project$class_clean[grep('oo',project$coarse_annotation_archr,ignore.case = T)]
project$class_clean[grep('oo',project$coarse_annotation_archr_h,ignore.case = T)]

project$coarse_annotation_archr[grep('oo',project$coarse_annotation_archr,ignore.case = T)]

# call peaks ###########
project = addGroupCoverages(ArchRProj = project, groupBy = "coarse_annotation_archr_clean")

geneAnnotation = readRDS('geneAnnotation_2020A.rds')

# run in /nfs/cellgeni/singularity/images/archr1.0.3_seurat4.2.1.sif
# as macs2 is not uninstallable
pathToMacs2 =  ArchR::findMacs2() # '/software/conda/users/pm19/macs2/bin/macs2'
project <- addReproduciblePeakSet(
  ArchRProj = project, 
  groupBy = "coarse_annotation_archr_clean", 
  pathToMacs2 = pathToMacs2,
  geneAnnotation = geneAnnotation
)

saveArchRProject(project,outputDirectory = 'work/archr2/02_annotate_out')
project = addPeakMatrix(project)
saveArchRProject(project,outputDirectory = 'work/archr2/02_annotate_out')

# peaks = ArchR::getPeakSet(project)
# saveRDS(peaks,'work/archr2/02_annotate_out/peaks.rds')
# write.csv(peaks,'work/archr2/02_annotate_out/peaks.csv')
# mtx = getMatrixFromProject(project,'PeakMatrix')
# saveRDS(mtx,'work/archr2/02_annotate_out/peaks2cell.rds')
mtx = readRDS('work/archr2/02_annotate_out/peaks2cell.rds')
for(cl in unique(mtx$class_clean)){
  print(cl)
  t = mtx[,mtx$class_clean==cl]
  saveRDS(t,paste0('work/archr2/03_per_class_annotation/',cl,'/peaks2cell.rds'))
}

# _reduce  dim on peaks ################
ArchR::getAvailableMatrices(project)


project <- addIterativeLSI(
  ArchRProj = project,
  useMatrix = "PeakMatrix", 
  name = "PeakMatrix_LSI", 
  iterations = 4, 
  clusterParams = list(resolution = c(2), sampleCells = 10000, maxClusters = 6, n.start = 10), 
  varFeatures = 25000, 
  dimsToUse = 1:30, 
  LSIMethod = 2
)


project = addHarmony(
  ArchRProj = project,
  reducedDims = "PeakMatrix_LSI",
  name = "PeakMatrix_LSI_harmony",
  groupBy = "Donor"
)

for(dim in c('PeakMatrix_LSI','PeakMatrix_LSI_harmony')){
  print(dim)
  for(r in c(0.2,1,2,10)){
    print(r)
    project <- addClusters(
      input = project,
      reducedDims = dim,
      method = "Seurat",
      name = paste0(dim,'_cl_',r),
      resolution = r, 
      #maxClusters = 25, 
      knnAssign = 20, 
      force = TRUE
    )
  }
  
  project <- addUMAP(
    ArchRProj = project,
    reducedDims = dim,
    name = paste0(dim,'_UMAP'),
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine"
  )
}
# doubl = filterDoublets(project)
# project$is_doublet = !(project$cellNames %in% doubl$cellNames)
# saveArchRProject(project,outputDirectory = 'work/archr2/02_annotate_out')


names(project@embeddings)
u = 'PeakMatrix_LSI_UMAP'
umap = project@embeddings[[u]]$df

par(mar=c(0,0,1,6))
visutils::plotVisium(umap,project$Donor,t='xy',cex=0.2,plot.legend = T,label.clusters = project$coarse_annotation_archr_clean)
visutils::plotVisium(umap,project$class_clean,t='xy',cex=0.2,plot.legend = T,label.clusters = T)
visutils::plotVisium(umap,project$is_doublet,t='xy',cex=0.2,plot.legend = T,label.clusters = T)
visutils::plotVisium(umap,project$coarse_annotation_archr_clean,t='xy',cex=0.2,plot.legend = T,label.clusters = T)
visutils::plotVisium(umap,project$DoubletScore,t='xy',cex=0.2,plot.legend = T,label.clusters = project$coarse_annotation_archr_clean)
visutils::plotVisium(umap,project$DoubletEnrichment,t='xy',cex=0.2,plot.legend = T,label.clusters = project$coarse_annotation_archr_clean,order.points.by.z = T,zfun = log1p)
visutils::plotVisium(umap,project$TSSEnrichment,t='xy',cex=0.2,plot.legend = T,label.clusters = project$coarse_annotation_archr_clean)
visutils::plotVisium(umap,project$nFrags,t='xy',cex=0.2,plot.legend = T,label.clusters = project$coarse_annotation_archr_clean)

u = 'PeakMatrix_LSI_harmony_UMAP'
umap = project@embeddings[[u]]$df
visutils::plotVisium(umap,project$is_doublet,t='xy',cex=0.2,plot.legend = T,label.clusters = T)
visutils::plotVisium(umap,project$Donor,t='xy',cex=0.2,plot.legend = T,label.clusters = paste(project$Donor,project$class_h))
visutils::plotVisium(umap,project$class_clean,t='xy',cex=0.2,plot.legend = T,label.clusters = T)

range(project$nFrags)

# check peak matrix #############
library(spam)
library(spam64)
project = loadArchRProject('work/archr/02_annotate_out')
mtx = getMatrixFromProject(project,'PeakMatrix')
# pc_mtx = schard::h5ad2list('work/pycistopic/results_pycistopic_call_peaks/combined.h5ad',use_spam = TRUE)
# dim(pc_mtx$X)
# 
# print(object.size(pc_mtx),unit='Gb')
peak = 'chr19:2248950-2249450'
# peak_inx = which(pc_mtx$var[,1]==peak)
# peak_cnts = pc_mtx$X[peak_inx,]
# peak_cnts = as.matrix(peak_cnts)[1,]
# names(peak_cnts) = paste0(pc_mtx$obs$sample_id,'#',pc_mtx$obs$barcode)
#saveRDS(peak_cnts,'work/pycistopic/results_pycistopic_call_peaks/peak_chr19_2248950-2249450_cov.rds')
peak_cnts = readRDS('work/pycistopic/results_pycistopic_call_peaks/peak_chr19_2248950-2249450_cov.rds')

i=findOverlaps(GRanges('chr19',IRanges(2248950 ,2249450)),rowRanges(mtx),type = 'any')@to
cmn = intersect(colnames(mtx),names(peak_cnts))
arch_cnts = assay(mtx,'PeakMatrix')[i,cmn]
# these are almost same
table(pycnt=peak_cnts[cmn],arch_cnts)


t = colSums(assay(mtx,'PeakMatrix'))[project$cellNames]
plot(project$ReadsInPeaks,t,pch='.')
hist(log2(t/project$nFrags/2),100)
plot(project$ReadsInPeaks,project$nFrags,pch='.')
hist((project$ReadsInPeaks/t),100)



# export #############
project = ArchR::loadArchRProject('work/archr2/02_annotate_out')
cols = c("Sample","TSSEnrichment","ReadsInTSS","ReadsInPromoter","ReadsInBlacklist","PromoterRatio","PassQC","NucleosomeRatio","nMultiFrags","nMonoFrags","nFrags","nDiFrags",
         "DoubletScore","DoubletEnrichment","BlacklistRatio","Age","Donor","TileMatrix_LSI_harmony_cl_0.2",
         "coarse_annotation_archr","class_clean","coarse_annotation_archr_clean","granulosa_refined_annotation",
         "ReadsInPeaks","FRIP","is_doublet","TileMatrix_LSI_cl_1","TileMatrix_LSI_cl_2","TileMatrix_LSI_cl_4","TileMatrix_LSI_cl_6","TileMatrix_LSI_cl_8","TileMatrix_LSI_cl_10")

gra = ArchR::loadArchRProject('work/archr2/04_per_class_clean/Granulosa')
project@cellColData[,'granulosa_refined_annotation'] =  NA
project@cellColData[gra$cellNames,'granulosa_refined_annotation'] = gra$coarse_annotation_final
table(project$class_clean,project$granulosa_refined_annotation)

umap = project@embeddings$TileMatrix_LSI_harmony_UMAP$df
umap = project@embeddings$TileMatrix_LSI_UMAP$df

TileMatrix_LSI_UMAP

plotVisium(umap,project$class_clean,t='xy',label.clusters = T,cex=0.5)
plotVisium(umap,project$Donor,t='xy',label.clusters = T,cex=0.5)

dir = paste0(project@projectMetadata$outputDirectory,'/export')
dir.create(dir)
write.csv(project@cellColData[,cols],paste0(dir,'/obs.csv'))
write.csv(project@embeddings$TileMatrix_LSI_UMAP$df,paste0(dir,'/TileMatrix_LSI_UMAP.csv'))
write.csv(project@embeddings$TileMatrix_LSI_harmony_UMAP$df,paste0(dir,'/TileMatrix_LSI_harmony_UMAP.csv'))
write.csv(project@reducedDims$TileMatrix_LSI$matSVD,paste0(dir,'/TileMatrix_LSI.csv'))
write.csv(project@reducedDims$TileMatrix_LSI_harmony$matDR,paste0(dir,'/TileMatrix_LSI_harmony.csv'))
getAvailableMatrices(project)
gsm = getMatrixFromProject(project,'GeneScoreMatrix')
saveRDS(gsm,'work/archr2/02_annotate_out/GeneScoreMatrix.rds')
#gsm = readRDS('work/archr2/02_annotate_out/GeneScoreMatrix.rds')
gsm_ = unwrapAssay(gsm)

Matrix::writeMM(gsm_,paste0(dir,'/genescores.mtx'))
write.csv(as.data.frame(rowData(gsm)),paste0(dir,'/genes.csv'))
write.csv(as.data.frame(colData(gsm))[,cols],paste0(dir,'/genescores_obs.csv'))


pm = getMatrixFromProject(project,'PeakMatrix')
pm_ = unwrapAssay(pm)
Matrix::writeMM(pm_,paste0(dir,'/PeakMatrix.mtx'))

write.csv(as.data.frame(rowRanges(pm)),paste0(dir,'/peaks.csv'))
write.csv(as.data.frame(colData(pm))[,cols],paste0(dir,'/PeakMatrix_obs.csv'))
