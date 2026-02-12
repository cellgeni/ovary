library(ArchR)
library(Seurat)
library(visutils)
source('actions/ovary/bin/archr_utils.R')

addArchRGenome("hg38")

addArchRThreads(threads = 2)
#project = loadArchRProject('work/archr_with_fetal/02_annotate')

# samples = read.csv('actions/samples.csv')
# samples$arrows = normalizePath(paste0('work/archr/01_arrows/arrows/',samples$sample_id,'/',samples$sample_id,'.arrow'))
# samples_f = read.csv('actions/samples_fetal.csv')
# samples_f$arrows = normalizePath(paste0('work/archr/01_arrows_fetal/arrows/',samples_f$sample_id,'/',samples_f$sample_id,'.arrow'))
# colnames(samples_f) = colnames(samples)
# samples = rbind(samples,samples_f)
# 
# all(file.exists(samples$arrows))
# 
# outdir = 'work/archr_with_fetal/02_annotate'
# 
# # annotate cells based on Valentinas approach
# # https://github.com/ventolab/Human-ReproductiveTract-Development-Atlas/blob/main/preprocessing/scatacseq/SingleSample_Analysis_ATAC_ArchR.ipynb
# 
# project <- ArchRProject(
#   ArrowFiles = samples$arrows,
#   outputDirectory = outdir,
#   copyArrows = TRUE
# )

project

geneAnnotation = myCreateGeneAnnotation('hg38','/software/cellgen/cellgeni/refdata_10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf')

agenes = ArchR::getGenes()
table(agenes$symbol %in% geneAnnotation$genes$symbol)
setdiff(agenes$symbol , geneAnnotation$genes$symbol)

geneAnnotation$genes = geneAnnotation$genes[geneAnnotation$genes$gene_type == 'protein_coding' | geneAnnotation$genes$symbol %in%  ArchR::getGenes()$symbol, ]
length(geneAnnotation$genes)
table(geneAnnotation$genes$gene_type)

project = addGeneScoreMatrix(project,genes = geneAnnotation$genes,force=TRUE)

# add metadata ######
samples_info = read.csv('data/Table1-dataset_metadata - scATACseq_SangerPediatric.csv',row.names = 1)[,c('Donor','Age')]
samples_info$sex = 'female'
samples_info$dataset = 'postnatal'
samples_info$technology = 'atac'

fetal_cells = read.csv('data/fetal/atac_metadata_with_multiome.csv')
fetal_cells$barcode = splitSub(fetal_cells$X,'-',2:3)
rownames(fetal_cells) = paste0(fetal_cells$sample,'#',fetal_cells$barcode)

fsamples =unique(fetal_cells[,c('individual','sample','sex','stage')])
rownames(fsamples) = fsamples$sample
fsamples = fsamples[rownames(fsamples) %in% project$Sample,]
fsamples$dataset = 'fetal'
fsamples$technology = ifelse(grepl('^HD',rownames(fsamples)),'atac','multiome')
fsamples$sample = NULL

colnames(fsamples)[1:3] = c('Donor','sex','Age')
samples_info = rbind(samples_info,fsamples)
all(project$Sample %in% rownames(samples_info))

project$Age = samples_info[project$Sample,'Age']
project$Donor = samples_info[project$Sample,'Donor']
project$sex = samples_info[project$Sample,'sex']
project$dataset = samples_info[project$Sample,'dataset']
project$technology = samples_info[project$Sample,'technology']

project$cell_type_fetal = fetal_cells[project$cellNames,'cell_type']
project$multiome_celltype_fetal = fetal_cells[project$cellNames,'multiome_celltype']

table(project$cell_type_fetal,paste(project$dataset,project$technology),useNA='always')
table(project$multiome_celltype_fetal,paste(project$dataset,project$technology),useNA='always')

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
# If you see downstream that you have subtle batch effects, another option is to add more LSI iterations and to start from a 
# lower intial clustering resolution as shown below. Additionally the number of variable features can be lowered to increase 
# focus on the more variable features. (https://www.archrproject.com/bookdown/iterative-latent-semantic-indexing-lsi.html)


project <- addIterativeLSI(
  ArchRProj = project,
  useMatrix = "TileMatrix", 
  name = "TileMatrix_LSI", 
  iterations = 4, 
  clusterParams = list(resolution = c(0.1,0.2,0.5), sampleCells = 50000, maxClusters = 6, n.start = 10), 
  varFeatures = 25000, 
  dimsToUse = 1:30, 
  LSIMethod = 2,
  force=TRUE
)


project <- addIterativeLSI(
  ArchRProj = project,
  useMatrix = "TileMatrix", 
  name = "TileMatrix_LSI_5000", 
  iterations = 4, 
  clusterParams = list(resolution = c(2), sampleCells = 10000, maxClusters = 6, n.start = 10), 
  varFeatures = 5000, 
  dimsToUse = 1:30, 
  LSIMethod = 2,
  force=TRUE
)


project <- addIterativeLSI(
  ArchRProj = project,
  useMatrix = "GeneScoreMatrix", 
  name = "GeneScoreMatrix_LSI", 
  iterations = 4, 
  clusterParams = list(resolution = c(2), sampleCells = 10000, maxClusters = 6, n.start = 10), 
  varFeatures = 5000, 
  dimsToUse = 1:30, 
  LSIMethod = 2,
  binarize = F,firstSelection = 'var',force = T
)


project = addHarmony(
  ArchRProj = project,
  reducedDims = "TileMatrix_LSI",
  name = "TileMatrix_LSI_harmony",
  groupBy = "Donor",
  force=TRUE
)

project = addHarmony(
  ArchRProj = project,
  reducedDims = "TileMatrix_LSI_5000",
  name = "TileMatrix_LSI_5000_harmony",
  groupBy = "Donor",
  force=TRUE
)


project = addHarmony(
  ArchRProj = project,
  reducedDims = "GeneScoreMatrix_LSI",
  name = "GeneScoreMatrix_LSI_harmony",
  groupBy = "Donor",
  force=TRUE
)




project$batch = paste0(project$dataset,'_',project$technology)

project = addHarmony(
  ArchRProj = project,
  reducedDims = "TileMatrix_LSI",
  name = "TileMatrix_LSI_harmony_batch",
  groupBy = "batch",
  force=TRUE
)

project = addHarmony(
  ArchRProj = project,
  reducedDims = "TileMatrix_LSI_5000",
  name = "TileMatrix_LSI_5000_harmony_batch",
  groupBy = "batch",
  force=TRUE
)



project = addHarmony(
  ArchRProj = project,
  reducedDims = "GeneScoreMatrix_LSI",
  name = "GeneScoreMatrix_LSI_harmony_batch",
  groupBy = "batch",
  force=TRUE
)



for(dim in c('TileMatrix_LSI_harmony_batch')){#names(project@reducedDims)){
  print(dim)
  for(r in c(4,6,10)){
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
  #   metric = "cosine",
  #   force=TRUE
  # )
}

project <- addModuleScore(project,
                            useMatrix = "GeneScoreMatrix",
                            name = "MS",
                            features = list(OSE=c('UPK3B','LRRN4','GATA4','LHX2')))

#pproj = loadArchRProject('work/archr2/02_annotate_out/')
project$class_clean_postnatal = NA
project@cellColData[pproj$cellNames,'class_clean_postnatal'] = pproj$class_clean

#saveArchRProject(project,outputDirectory = 'work/archr_with_fetal/02_annotate')

# gs = ArchR::getMatrixFromProject(project,useMatrix = 'GeneScoreMatrix')
# saveRDS(gs,'work/archr_with_fetal/02_annotate/GeneScoreMatrix.rds')
project = loadArchRProject('work/archr_with_fetal/02_annotate')
gs = readRDS('work/archr_with_fetal/02_annotate/GeneScoreMatrix.rds')
mtx = assay(gs,'GeneScoreMatrix')[,project$cellNames]
rownames(mtx) = rowData(gs)$name
colnames(mtx) = project$cellNames

c('UPK3B','LRRN4','GATA4','LHX2') %in% rowData(gs)$name

# markerGenes = c('LHX2')
# p <- plotEmbedding(
#   ArchRProj = project, 
#   colorBy = "GeneScoreMatrix", 
#   name = markerGenes, 
#   embedding = "TileMatrix_LSI_UMAP",
#   quantCut = c(0.01, 0.95),
#   imputeWeights = NULL
# )
# p

lrrn4 = mtx['LRRN4',]
lhx2 = mtx['LHX2',]
upk3b = mtx['UPK3B',]
gata2 = mtx['GATA2',]
gata4 = mtx['GATA4',]
l = project$class_clean_postnatal
l[is.na(l)] = paste0("F",ifelse(project$technology[is.na(l)]=='atac','a','m'),':',project$cell_type_fetal[is.na(l)])
f = !startsWith(l,'F')
l[f] = paste0('P:',l[f])
ose = project$cell_type_fetal == 'CoelEpi'
ose[is.na(ose)] = F
table(ose)
f=function(x,v='-'){ifelse(is.na(x),v,x)}


table(project$cell_type_fetal,project$batch)


png('work/archr_with_fetal/02_annotate/figures/umaps2.png',w=9*7,h=7*4.5,unit='in',res = 200)
par(mfcol=c(7,9),mar=c(1,1,1,15),oma=c(0,0,1,0))
for(u in names(project@embeddings)){
  umap = project@embeddings[[u]]$df
  plotVisium(umap,project$batch,t='xy',cex=0.2,randomize.points = T,main=u)
  plotVisium(umap,l,t='xy',cex=0.2,label.clusters = l,order.points.by.z = T,main=u)
  plotVisium(umap,ose,t='xy',cex=ose+0.2,label.clusters = l,order.points.by.z = T,main=u,legend.args = list(title='is CoelEpi'))
  plotVisium(umap,trimQ(lhx2,0.999),t='xy',cex=0.2,label.clusters = l,order.points.by.z = T,main=u,legend.args = list(title='LHX2'),zfun=log1p)
  plotVisium(umap,trimQ(upk3b,0.999),t='xy',cex=0.2,label.clusters = l,order.points.by.z = T,main=u,legend.args = list(title='UPK3B'),zfun=log1p)
  plotVisium(umap,trimQ(gata2,0.999),t='xy',cex=0.2,label.clusters = l,order.points.by.z = T,main=u,legend.args = list(title='GATA2'),zfun=log1p)
  plotVisium(umap,trimQ(project$MS.OSE,0.995),t='xy',cex=0.2,label.clusters = l,randomize.points = T,main=u,legend.args = list(title='OSE score'))
  
  plotVisium(umap,l,t='xy',cex=0.2,label.clusters = l,order.points.by.z = T,spot.filter = project$batch!='postnatal_atac',main=u)
}
dev.off()

u='TileMatrix_LSI_harmony_batch_UMAP'
umap = project@embeddings[[u]]$df
par(mfcol=c(1,1),mar=c(1,1,1,10),oma=c(0,0,1,0))
plotVisium(umap,l,t='xy',cex=0.2,label.clusters = l,order.points.by.z = T,
           spot.filter = project$batch!='postnatal_atac' & project$cell_type_fetal %in% c('CoelEpi','preGranulosa'),
           main=u)

plotVisium(umap,project$TileMatrix_LSI_harmony_batch_cl_6,t='xy',cex=0.2,label.clusters = T,order.points.by.z = T,main=u)

table(project$TileMatrix_LSI_harmony_batch_cl_6,project$TileMatrix_LSI_harmony_batch_cl_10)

ff = function(v,f){
  r = split(v,f)
  r[order(sapply(r,mean),decreasing = T)]
}

table(project$TileMatrix_LSI_harmony_batch_cl_6,l)[,c('Fa:CoelEpi','Fm:CoelEpi')]
t(table(project$TileMatrix_LSI_harmony_batch_cl_6,l)[c('C20','C13','C14'),])
# C13 Fa:CoelEpi +
# C20 Fm:CoelEpi +
# C14 P:Granulosa (looks like OSE)

par(mar=c(15,4,1,1))
boxplot(ff(project$MS.OSE ,project$TileMatrix_LSI_harmony_batch_cl_6),las=2)
boxplot(ff(lhx2 ,project$TileMatrix_LSI_harmony_batch_cl_6),las=2)
boxplot(ff(upk3b ,project$TileMatrix_LSI_harmony_batch_cl_6),las=2)
boxplot(ff(lrrn4 ,project$TileMatrix_LSI_harmony_batch_cl_6),las=2)
boxplot(ff(gata4 ,project$TileMatrix_LSI_harmony_batch_cl_6),las=2)
project$potential_OSE = project$TileMatrix_LSI_harmony_batch_cl_6=='C14'
table(l,project$potential_OSE)
#saveArchRProject(project,outputDirectory = 'work/archr_with_fetal/02_annotate')
