library(ArchR)
library(Seurat)
library(visutils)
source('actions/ovary/bin/archr_utils.R')

addArchRGenome("hg38")

addArchRThreads(threads = 1)
#project = loadArchRProject('work/archr_with_fetal_atac/02_annotate')

samples = read.csv('actions/samples.csv')
samples$arrows = normalizePath(paste0('work/archr/01_arrows/arrows/',samples$sample_id,'/',samples$sample_id,'.arrow'))
samples_f = read.csv('actions/samples_fetal.csv')
samples_f$arrows = normalizePath(paste0('work/archr/01_arrows_fetal/arrows/',samples_f$sample_id,'/',samples_f$sample_id,'.arrow'))
colnames(samples_f) = colnames(samples)
samples = rbind(samples,samples_f)
samples$technology = ifelse(grepl('^FCA',samples$sample_id),'multiome','atac')
# only atac
samples = samples[samples$technology!='multiome',]

all(file.exists(samples$arrows))

outdir = 'work/archr_with_fetal_atac/02_annotate'

# annotate cells based on Valentinas approach
# https://github.com/ventolab/Human-ReproductiveTract-Development-Atlas/blob/main/preprocessing/scatacseq/SingleSample_Analysis_ATAC_ArchR.ipynb

project <- ArchRProject(
  ArrowFiles = samples$arrows,
  outputDirectory = outdir,
  copyArrows = TRUE
)

project

geneAnnotation = myCreateGeneAnnotation('hg38','/software/cellgen/cellgeni/refdata_10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf')

# with links matrix it too large
geneAnnotation$genes = geneAnnotation$genes[geneAnnotation$genes$gene_type == 'protein_coding' | geneAnnotation$genes$symbol %in%  ArchR::getGenes()$symbol, ]

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
project$batch = paste0(project$dataset,'_',project$technology)


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
  name = "Tile_LSI", 
  iterations = 4, 
  clusterParams = list(resolution = c(0.1,0.2,0.5), sampleCells = 50000, maxClusters = 6, n.start = 10), 
  varFeatures = 25000, 
  dimsToUse = 2:40, 
  LSIMethod = 2,
  force=TRUE
)
reducedDims = "Tile_LSI"
for(batch in c("Donor","batch",'Sample')){
  project = addHarmony(
    ArchRProj = project,
    reducedDims = reducedDims,
    name = paste0(reducedDims,'_',batch),
    groupBy = batch,
    force=TRUE
  )
}


for(dim in 'Tile_LSI_batch'){#names(project@reducedDims)){
  print(dim)
  for(r in c(4,10)){#c(0.2,1,2)){
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

#saveArchRProject(project,outputDirectory = 'work/archr_with_fetal_atac/02_annotate')

#gs = ArchR::getMatrixFromProject(project,useMatrix = 'GeneScoreMatrix')
#saveRDS(gs,'work/archr_with_fetal_atac/02_annotate/GeneScoreMatrix.rds')

gs = readRDS('work/archr_with_fetal_atac/02_annotate/GeneScoreMatrix.rds')
mtx = assay(gs,'GeneScoreMatrix')[,project$cellNames]
rownames(mtx) = rowData(gs)$name
colnames(mtx) = project$cellNames


lhx2 = mtx['LHX2',]
upk3b = mtx['UPK3B',]
gata2 = mtx['GATA2',]
l = project$class_clean_postnatal
l[is.na(l)] = paste0("F:",project$cell_type_fetal[is.na(l)])
f = !startsWith(l,'F')
l[f] = paste0('P:',l[f])
ose = project$cell_type_fetal == 'CoelEpi'
ose[is.na(ose)] = F
table(ose)
f=function(x,v='-'){ifelse(is.na(x),v,x)}


table(l)

table(project$cell_type_fetal,project$batch)
ncol = length(project@embeddings)

png('work/archr_with_fetal_atac/02_annotate/figures/umaps2.png',h=ncol*4.5,w=7*8,unit='in',res = 200)
par(mfrow=c(ncol,8),mar=c(1,1,1,15),oma=c(0,0,1,0))
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

u='Tile_LSI_batch_UMAP'
umap = project@embeddings[[u]]$df
plotVisium(umap,project$Tile_LSI_batch_cl_10,t='xy',cex=0.2,label.clusters = T,order.points.by.z = T,main=u)
table(project$Tile_LSI_batch_cl_10,l)[,'F:CoelEpi']
table(l,project$Tile_LSI_batch_cl_10=='C8')


ff = function(v,f){
  r = split(v,f)
  r[order(sapply(r,mean),decreasing = T)]
}

par(mar=c(15,4,1,1))
boxplot(ff(project$MS.OSE ,l),las=2)
boxplot(ff(project$MS.OSE ,project$Tile_LSI_batch_cl_10),las=2)
boxplot(ff(lhx2 ,project$Tile_LSI_batch_cl_10),las=2)
boxplot(ff(upk3b ,project$Tile_LSI_batch_cl_10),las=2)


# integrate ############
rna_cells = read.csv('data/fetal/scRNAseq_directannotation_perlineagescvi_withFetal.csv',row.names = 1)
z=splitSub(rownames(rna_cells),'-',1)
rna_cells$sample_id = substr(z,1,nchar(z)-17)
rna_cells$barcode = substr(z,nchar(z)-15,10000)

rna = schard::h5ad2seurat('data/ref_rawcounts.h5ad')
rna$coarse_annotation

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
visutils::plotVisium(umap,project$class_clean,t='xy',cex=0.1,label.clusters = T,randomize.points = T)
visutils::plotVisium(umap,project$class_h=='PV',t='xy',order.points.by.z = T,cex=(project$class_h=='PV')+0.1)
visutils::plotVisium(umap,project$TileMatrix_LSI_cl_1=='C3',t='xy',order.points.by.z = T,cex=(project$TileMatrix_LSI_cl_1=='C3')+0.1)

# annotated per coarse celltypes groups #########
project = ArchR::loadArchRProject('work/archr2/02_annotate_out')
path = 'work/archr2/03_per_class_annotation'
rna = schard::h5ad2seurat('data/ref_rawcounts.h5ad')
rna$class = visutils::splitSub(rna$coarse_annotation,'_',1)
rna$class[rna$class %in% c('Theca','Stroma','FibC3')] = 'Mesenchymal'

#dir.create(path)
setdiff(project$class_clean,rna$class)
ArchR::addArchRVerbose(FALSE)

for(class  in unique(project$class_clean)){
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
  
  # _integrate ############
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
project = ArchR::loadArchRProject('work/archr_with_fetal_atac/02_annotate')
cols = colnames(project@cellColData)


dir = paste0(project@projectMetadata$outputDirectory,'/export')
dir.create(dir)
write.csv(project@cellColData[,cols],paste0(dir,'/obs.csv'))
for(u in names(project@embeddings))
  write.csv(project@embeddings[[u]]$df,paste0(dir,'/',u,'.csv'))
for(d in names(project@reducedDims))
  write.csv(project@reducedDims[[d]]$matDR,paste0(dir,'/',d,'.csv'))

gsm = getMatrixFromProject(project,'GeneScoreMatrix')
saveRDS(gsm,paste0(project@projectMetadata$outputDirectory,'/GeneScoreMatrix.rds'))
#gsm = readRDS(paste0(project@projectMetadata$outputDirectory,'/GeneScoreMatrix.rds'))
gsm_ = unwrapAssay(gsm)

Matrix::writeMM(gsm_,paste0(dir,'/genescores.mtx'))
write.csv(as.data.frame(rowData(gsm)),paste0(dir,'/genes.csv'))
write.csv(as.data.frame(colData(gsm))[,cols],paste0(dir,'/genescores_obs.csv'))


prj = ArchR::loadArchRProject('work/archr2/02_annotate_out/')
peaks = ArchR::getPeakSet(prj)
project = addPeakSet(project,peaks)
project = addPeakMatrix(project)
saveArchRProject(project,outputDirectory = 'work/archr_with_fetal_atac/02_annotate')

getAvailableMatrices(project)
pm = getMatrixFromProject(project,'PeakMatrix')
pm_ = unwrapAssay(pm)
Matrix::writeMM(pm_,paste0(dir,'/PeakMatrix.mtx'))

write.csv(as.data.frame(rowRanges(pm)),paste0(dir,'/peaks.csv'))
write.csv(as.data.frame(colData(pm))[,cols],paste0(dir,'/PeakMatrix_obs.csv'))


