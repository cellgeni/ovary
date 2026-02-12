# This script is deprictated, see 01_archr_qc.R
library(ArchR)
library(Seurat)
library(visutils)

addArchRGenome("hg38")
addArchRThreads(threads = 4)

samples = read.csv('actions/samples.csv')
samples$arrows = normalizePath(paste0('work/archr/01_arrows/arrows/',samples$sample_id,'/',samples$sample_id,'.arrow'))
file.exists(samples$arrows)

#samples = samples[1,]
outdir = 'work/archr/02_annotate'

# annotate cells based on Valentinas approach
# https://github.com/ventolab/Human-ReproductiveTract-Development-Atlas/blob/main/preprocessing/scatacseq/SingleSample_Analysis_ATAC_ArchR.ipynb

project <- ArchRProject(
  ArrowFiles = samples$arrows,
  outputDirectory = outdir,
  copyArrows = TRUE
)

project

# add external annotation ###########
scanvi_ann = schard::h5ad2data.frame('work/scanvi_out/combined_gene_matrix_plus_ref_clean_scanvi_raw.h5ad','/obs')
scanvi_ann = scanvi_ann[scanvi_ann$dataset=='atac',]
rownames(scanvi_ann) = paste0(scanvi_ann$sample,'#',visutils::splitSub(rownames(scanvi_ann),'[:_]',5,fixed = F))
table(project$cellNames %in% rownames(scanvi_ann))
table(rownames(scanvi_ann) %in% project$cellNames)

cmn = intersect(rownames(scanvi_ann) , project$cellNames)
project@cellColData[,'coarse_annotation_bbknn'] = NA
project@cellColData[cmn,'coarse_annotation_bbknn'] = scanvi_ann[cmn,'coarse_annotation_bbknn']

samples_info = read.csv('data/Table1-dataset_metadata - scATACseq_SangerPediatric.csv',row.names = 1)
project$Age = samples_info[project$Sample,'Age']
project$Donor = samples_info[project$Sample,'Donor']


# add peaks from pycistopic ############
peaks = read.table('work/pycistopic/results_pycistopic_call_peaks/consensus_peaks.bed')
peaks = GRanges(peaks$V1,ranges = IRanges(start=peaks$V2,end = peaks$V3))
project = addPeakSet(project,peaks)
project = addPeakMatrix(project)
getAvailableMatrices(project)

# m = getMatrixFromProject(project,'PeakMatrix')
# plot(project$ReadsInPeaks,colSums(assay(m,'PeakMatrix')),pch='.')
# hist((project$ReadsInPeaks/colSums(assay(m,'PeakMatrix'))))

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
  name = "whole_tile_LSI", 
  iterations = 2, 
  clusterParams = list(resolution = c(2), sampleCells = 10000, maxClusters = 6, n.start
                       = 10), 
  varFeatures = 25000, 
  dimsToUse = 1:30, 
  LSIMethod = 2
)


project <- addIterativeLSI(
  ArchRProj = project,
  useMatrix = "PeakMatrix", 
  name = "whole_peak_LSI", 
  iterations = 2, 
  clusterParams = list(resolution = c(2), sampleCells = 10000, maxClusters = 6, n.start
                       = 10), 
  varFeatures = 25000, 
  dimsToUse = 1:30, 
  LSIMethod = 2
)


project = addHarmony(
  ArchRProj = project,
  reducedDims = "whole_peak_LSI",
  name = "whole_peak_LSI_harmony",
  groupBy = "Donor"
)

for(dim in names(project@reducedDims)){
  print(dim)
  project <- addClusters(
    input = project,
    reducedDims = dim,
    method = "Seurat",
    name = paste0(dim,'_cl_1'),
    resolution = 1, 
    #maxClusters = 25, 
    knnAssign = 20, 
    force = TRUE
  )
  
  
  project <- addUMAP(
    ArchRProj = project, 
    reducedDims = dim, 
    name = paste0(dim,'_UMAP'), 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
  )
}

# integrate ############
rna = schard::h5ad2seurat('data/ref_rawcounts.h5ad')
rna$coarse_annotation

project <- addGeneIntegrationMatrix(
  ArchRProj = project, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "whole_peak_LSI_harmony",
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
  reducedDims = "whole_peak_LSI",
  seRNA = rna,
  addToArrow = FALSE,
  force = TRUE,
  groupRNA = 'coarse_annotation',
  nameCell = "coarse_annotation_archr_cell",
  nameGroup = "coarse_annotation_archr",
  nameScore = "coarse_annotation_archr_score"
)
saveArchRProject(project,outputDirectory = 'work/archr/02_annotate_out')

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

# export peak*cell matrix ###########
#obs = schard::h5ad2data.frame('data/ref_rawcounts.h5ad','obs')

# project = loadArchRProject('work/archr/02_annotate_out')
# mtx = getMatrixFromProject(project,'PeakMatrix')
# granuloasa = readRDS('work/archr/03_granulosa_clean_v2/Save-ArchR-Project.rds')
# gcells = granuloasa@cellColData
# mtx$coarse_annotation_granulosa_refined = NA
# mtx@colData[rownames(gcells),'coarse_annotation_granulosa_refined'] = gcells$coarse_annotation_archr_clean
# 
# table(mtx$coarse_annotation_bbknn,mtx$coarse_annotation_granulosa_refined)
# table(mtx$coarse_annotation_archr,mtx$coarse_annotation_granulosa_refined)
# table(splitSub(mtx$coarse_annotation_bbknn,'_',1),splitSub(mtx$coarse_annotation_archr_h,'_',1))
# 
# 
# saveRDS(mtx,'work/archr/02_annotate_out/PeakMatrix.rds')
mtx = readRDS('work/archr/02_annotate_out/PeakMatrix.rds')
mtx_m = t(assay(mtx,'PeakMatrix'))
obs = as.data.frame(mtx@colData)
var = as.data.frame(mtx@rowRanges)
all(rownames(obs) == rownames(mtx_m))

Matrix::writeMM(mtx_m,'work/archr/02_annotate_out/PeakMatrix.mtx')
write.csv(obs,'work/archr/02_annotate_out/PeakMatrix_obs.csv')
write.csv(var,'work/archr/02_annotate_out/PeakMatrix_var.csv')

# import scanpy as sc
# import anndata
# import pandas as pd
# 
# X = sc.read_mtx('PeakMatrix.mtx')
# obs = pd.read_csv('PeakMatrix_obs.csv',index_col=0)
# var = pd.read_csv('PeakMatrix_var.csv',index_col=0)
# var.index = var.seqnames+':'+var['start'].astype(str)+"-"+var['end'].astype(str)
# adata = anndata.AnnData(X=X.X,obs=obs,var=var)
# adata.write_h5ad('PeakMatrix.h5ad')
