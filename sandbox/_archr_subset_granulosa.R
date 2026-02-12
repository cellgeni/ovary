library(ArchR)
library(Seurat)
library(visutils)
library(SummarizedExperiment)
source('actions/ovary/bin/plotATACCoverage.R')
source('actions/ovary/bin/archr_utils.R')

addArchRGenome("hg38")
addArchRThreads(threads = 4)
ord = c('Granulosa_sq','Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml')
#granulosa = loadArchRProject('work/archr/03_granulosa')
gene.descr = readRDS('/nfs/cellgeni/pasham/code/nf-scsajr/ref/human_2020A/functional_annotation/gene.descr.rds')

gsm = getMatrixFromProject(granulosa,'GeneScoreMatrix')[,granulosa$cellNames]
gsm_genes = rowData(gsm)
gsm = assay(gsm,'GeneScoreMatrix')
rownames(gsm) = gsm_genes$name


# links peaks to genes ########
gex = schard::h5ad2list('data/ref_rawcounts.h5ad')

peaks = read.table('work/pycistopic/results_pycistopic_call_peaks/consensus_peaks.bed')
rownames(peaks) = paste0(peaks$V1,':',peaks$V2,'-',peaks$V3)
peaksg = GRanges(peaks$V1,ranges = IRanges(start=peaks$V2,end = peaks$V3))
tss = ifelse(gene.descr$strand==1,gene.descr$start,gene.descr$end)
genes_tss = GRanges(paste0('chr',gene.descr$chr),ranges = IRanges(start = tss,end = tss))
p2g_ = findOverlaps(peaksg,genes_tss,maxgap=250000,ignore.strand=TRUE)

p2g = data.frame(peak = rownames(peaks)[p2g_@from],gene_id=gene.descr$ens_id[p2g_@to])
p2g$gene_name = gene.descr[p2g$gene_id,'name']
p2g$dist2tss = gene.descr$strand[p2g_@to]*(tss[p2g_@to] - peaks$V2[p2g_@from]/2-peaks$V3[p2g_@from]/2)
hist(p2g$dist,100)
p2g[p2g$peak=='chr19:2248950-2249450',]

peaks_with_genes = do.call(rbind,lapply(split(p2g,p2g$peak),function(x){
  x = x[order(abs(x$dist2tss)),]
  data.frame(nearest_gene_id = x$gene_id[1],
             nearest_gene_name = x$gene_name[1],
             nearest_gene_dist2tss = x$dist2tss[1],
             all_gene_ids = paste(x$gene_id,collapse = ','),
             all_gene_names = paste(x$gene_name,collapse = ',')
             )
}))
peaks_with_genes_ = peaks_with_genes[rep(1,nrow(peaks)),]
rownames(peaks_with_genes_) = rownames(peaks)
peaks_with_genes_[,] = NA
peaks_with_genes_[rownames(peaks_with_genes),] = peaks_with_genes
peaks_with_genes = peaks_with_genes_
peaks_with_genes[1:2,]
peaks_with_genes['chr19:2248950-2249450',]
# write.csv(peaks_with_genes,'work/pycistopic/results_pycistopic_call_peaks/consensus_peaks_with_genes.csv')
# write.csv(p2g,'work/pycistopic/results_pycistopic_call_peaks/consensus_peaks2genes.csv')

# check arch emb/ann ##########
# looks like harmony embedding is bad - it splits granulosa cells into two clusters
project = loadArchRProject('work/archr/02_annotate_out')
umap=project@embeddings$whole_tile_LSI_UMAP$df
ct = project$coarse_annotation_archr

ctno = as.character(as.numeric(factor(ct)))
f = nchar(ctno) == 1
ctno[f] = paste0('0',ctno[f])
par(mfrow=c(1,2),mar=c(1,1,1,6))
plotVisium(umap,paste(ctno,ct),t='xy',cex=0.2,label.clusters = ctno,legend.args = list(ncol=2))
plotVisium(umap,grepl('Granulosa',ct),t='xy',cex=0.2,label.clusters = ctno)
plotVisium(umap,project$Sample,t='xy',cex=0.2)

par(mfrow=c(2,2),mar=c(1,1,1,1))
plotVisium(umap,grepl('Granulosa',project$coarse_annotation_archr),t='xy',cex=0.2,label.clusters = ctno)
plotVisium(umap,grepl('Granulosa',project$coarse_annotation_archr_h),t='xy',cex=0.2,label.clusters = ctno)
plotVisium(umap,grepl('Granulosa',project$coarse_annotation_bbknn),t='xy',cex=0.2,label.clusters = ctno)

table(scanvi=startsWith(project$coarse_annotation_bbknn,'Granulosa'),
      archr=startsWith(project$coarse_annotation_archr,'Granulosa'),useNA='always')

table(scanvi=startsWith(project$coarse_annotation_bbknn,'Granulosa'),
      archr=startsWith(project$coarse_annotation_archr_h,'Granulosa'),useNA='always')

table(archr_h=startsWith(project$coarse_annotation_archr_h,'Granulosa'),
      archr=startsWith(project$coarse_annotation_archr,'Granulosa'),useNA='always')


normMtx = function(m,fun=min){
  rows = rowSums(m)
  cols = colSums(m)
  for(r in seq_along(rows)){
    for(c in seq_along(cols)){
      m[r,c] = m[r,c]/fun(rows[r],cols[c]) 
    }
  }
  m
}
f = grepl('^Granulosa',project$coarse_annotation_bbknn) & grepl('^Granulosa',project$coarse_annotation_archr)
t = table(project$coarse_annotation_bbknn[f],project$coarse_annotation_archr[f])

# f = grepl('^Granulosa',project$coarse_annotation_bbknn) & grepl('^Granulosa',project$coarse_annotation_archr_h)
# t = table(project$coarse_annotation_bbknn[f],project$coarse_annotation_archr_h[f])
o = names(sort(rowSums(t),decreasing = T))
t = t[o,o[o %in% colnames(t)]]
par(mfrow=c(1,1),mar=c(10,10,1,1))

visutils::imageWithText(normMtx(t,max),t)


# subset granilosa #############
granulosa <- subsetArchRProject(
  ArchRProj = project,
  cells = project$cellNames[f],
  outputDirectory = "work/archr/03_granulosa",
  dropCells = TRUE,
  force = TRUE
)

# _check peak matrix ##################
getAvailableMatrices(granulosa)
mtx = getMatrixFromProject(granulosa,'PeakMatrix')
peak = 'chr19:2248950-2249450'
peak_cnts = readRDS('work/pycistopic/results_pycistopic_call_peaks/peak_chr19_2248950-2249450_cov.rds')
i=findOverlaps(GRanges('chr19',IRanges(2248950 ,2249450)),rowRanges(mtx),type = 'any')@to
cmn = intersect(colnames(mtx),names(peak_cnts))
arch_cnts = assay(mtx,'PeakMatrix')[i,cmn]

# complete bullshit!
table(pycnt=peak_cnts[cmn],arch_cnts)

# lets add it again!
peaks = read.table('work/pycistopic/results_pycistopic_call_peaks/consensus_peaks.bed')
peaks = GRanges(peaks$V1,ranges = IRanges(start=peaks$V2,end = peaks$V3))
granulosa = addPeakSet(granulosa,peaks,force = TRUE)
granulosa = addPeakMatrix(granulosa,force=TRUE)



# dim reduction #########
granulosa <- addIterativeLSI(
  ArchRProj = granulosa,
  useMatrix = "TileMatrix", 
  name = "g_tile_iLSI", 
  iterations = 2, 
  clusterParams = list(resolution = c(2), sampleCells = 10000, maxClusters = 6, n.start = 10), 
  varFeatures = 25000, 
  dimsToUse = 1:30, 
  LSIMethod = 2
)


granulosa <- addClusters(
  input = granulosa,
  reducedDims = "g_tile_iLSI",
  method = "Seurat",
  name = "g_tile_Clusters_1",
  resolution = 1, 
  knnAssign = 20, 
  force = TRUE
)

granulosa <- addUMAP(
  ArchRProj = granulosa, 
  reducedDims = "g_tile_iLSI", 
  name = "g_tile_UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)


# dim reduction on peak matrix #########
granulosa <- addIterativeLSI(
  ArchRProj = granulosa,
  useMatrix = "PeakMatrix", 
  name = "g_peak_iLSI", 
  iterations = 2, 
  clusterParams = list(resolution = c(2), sampleCells = 10000, maxClusters = 6, n.start = 10), 
  varFeatures = 25000, 
  dimsToUse = 1:30, 
  LSIMethod = 2
)

granulosa <- addUMAP(
  ArchRProj = granulosa, 
  reducedDims = "g_peak_iLSI", 
  name = "g_peak_UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

for(r in c(1,2,3,10,0.5)){
  granulosa <- addClusters(
    input = granulosa,
    reducedDims = "g_peak_iLSI",
    method = "Seurat",
    name = paste0("g_peak_Clusters_",r),
    resolution = r, 
    knnAssign = 20, 
    force = TRUE
  )
}


# harmony ##########
granulosa <- addHarmony(
  ArchRProj = granulosa,
  reducedDims = "g_peak_iLSI",
  name = "g_peak_iLSI_harmony",
  groupBy = "Donor"
)

granulosa <- addUMAP(
  ArchRProj = granulosa, 
  reducedDims = "g_peak_iLSI_harmony", 
  name = "g_peak_harmony_UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

granulosa <- addClusters(
  input = granulosa,
  reducedDims = "g_peak_iLSI_harmony",
  method = "Seurat",
  name = "g_peak_harmony_Clusters_10",
  resolution = 10, 
  knnAssign = 20, 
  force = TRUE
)

# umaps on while red dims ##
for(d in c('whole_tile_LSI','whole_peak_LSI','whole_peak_LSI_harmony')){
  granulosa <- addUMAP(
    ArchRProj = granulosa, 
    reducedDims = d, 
    name = paste0(d,"_UMAP"), 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
  )
}

# refine annotation ###########
cols = char2col(granulosa$coarse_annotation_bbknn)
pdf('actions/ovary/figures/arch_comp_umaps.pdf',w=6*5,h=4*3)
par(mfcol=c(4,6),mar=c(1,1,1,14))
for(e in names(granulosa@embeddings)){
  for(ann in c('coarse_annotation_bbknn','coarse_annotation_archr','Sample','nFrags')){
    cols_ = num2col
    if(startsWith(ann,'coarse'))
      cols_ = cols
    visutils::plotVisium(granulosa@embeddings[[e]]$df,cex=0.8,z2col=cols_,
                         as.data.frame(granulosa@cellColData)[,ann],t='xy',main=e,
                         legend.args = list(title=ann))
  }
}
dev.off()

# whole_peak_LSI looks best, lets use it ###
granulosa <- addClusters(
  input = granulosa,
  reducedDims = "g_peak_iLSI",
  method = "Seurat",
  name = paste0("g_peak_Clusters_",10),
  resolution = 10, 
  knnAssign = 20, 
  force = TRUE
)


pdf('actions/ovary/figures/archr_granulosa_ann.pdf',w=30,h=16)
par(mfrow=c(3,4),mar=c(1,1,2,10),oma=c(0,0,0,6))
umap = granulosa@embeddings$g_peak_UMAP$df
cex=0.9
cols = visutils::char2col(granulosa$coarse_annotation_bbknn)
#visutils::plotVisium(umap,granulosa$peak_Clusters_10,t='xy',label.clusters = T,xlab='',ylab='',main='Clusters',cex=cex)
visutils::plotVisium(umap,granulosa$g_peak_Clusters_10,t='xy',label.clusters = T,xlab='',ylab='',main='Clusters',cex=cex)
visutils::plotVisium(umap,granulosa$coarse_annotation_bbknn,t='xy',label.clusters = T,xlab='',ylab='',main='scANvi (pycistopic GS)',z2col = cols,show.cluster.sizes = T,randomize.points = T,cex=cex)
visutils::plotVisium(umap,granulosa$coarse_annotation_archr,t='xy',label.clusters = T,xlab='',ylab='',main='ArchR label transfer',z2col = cols,show.cluster.sizes = T,randomize.points = T,cex=cex)
visutils::plotVisium(umap,granulosa$Sample,t='xy',label.clusters = F,xlab='',ylab='',main='Sample',cex=cex)
visutils::plotVisium(umap,granulosa$Donor,t='xy',label.clusters = F,xlab='',ylab='',main='Donor',cex=cex)
visutils::plotVisium(umap,granulosa$DoubletScore,t='xy',label.clusters = F,xlab='',ylab='',main='DoubletScore',cex=cex,order.points.by.z = T)
visutils::plotVisium(umap,granulosa$TSSEnrichment,t='xy',label.clusters = F,xlab='',ylab='',main='TSSEnrichment',cex=cex,order.points.by.z = T)
visutils::plotVisium(umap,granulosa$nFrags,t='xy',label.clusters = F,xlab='',ylab='',main='nFrags',cex=cex,order.points.by.z = T)
visutils::plotVisium(umap,granulosa$PromoterRatio,t='xy',label.clusters = F,xlab='',ylab='',main='PromoterRatio',cex=cex,order.points.by.z = T)
visutils::plotVisium(umap,granulosa$Age,t='xy',label.clusters = F,xlab='',ylab='',main='Age',cex=cex,order.points.by.z = T)
dev.off()

table(granulosa$Sample,granulosa$coarse_annotation_bbknn)


annByMaj = function(cl,ann){
  t = table(cl,ann)
  inx = apply(t,1,which.max)
  r = data.frame(cl=rownames(t),ann=colnames(t)[inx])
  r$cnt = as.numeric(table(cl)[r$cl])
  r[order(r$cl),]
}

(tb = table(granulosa$coarse_annotation_bbknn,granulosa$g_peak_Clusters_10))
(ta = table(granulosa$coarse_annotation_archr,granulosa$g_peak_Clusters_10))
sort(ta['Granulosa_sq_transitioning',]/colSums(ta)*100)
sort(tb['Granulosa_sq_transitioning',]/colSums(tb)*100)
cl2ann = annByMaj(granulosa$g_peak_Clusters_10,granulosa$coarse_annotation_bbknn)
cl2ann = cbind(cl2ann,ann_archr=annByMaj(granulosa$g_peak_Clusters_10,granulosa$coarse_annotation_archr)$ann)
rownames(cl2ann) = cl2ann$cl
#cl2ann[c('C15','C18'),'ann_archr'] = 'Granulosa_sq_transitioning' # clusters with highest Granulosa_sq_transitioning fraction
cl2ann[c('C2','C3'),'ann_archr'] = 'Granulosa_sq_transitioning' # clusters with highest Granulosa_sq_transitioning fraction
cl2ann[c('C23'),'ann_archr'] = 'Granulosa_AMH_early'

# for now I do not use harmony as it seems to be overcorrected (let try run harmony on whole obj)
# and I'm using arch annotation here as it give bit more of AMH_early
granulosa$coarse_annotation_archr_clean = cl2ann[granulosa$g_peak_Clusters_10,'ann_archr']
umap = granulosa@embeddings$g_peak_UMAP$df
par(mar=c(1,1,2,10))
visutils::plotVisium(umap,granulosa$coarse_annotation_archr_clean,t='xy',label.clusters = T,xlab='',ylab='',main='scANvi (pycistopic GS)',show.cluster.sizes = T,cex=cex,order.points.by.z = T)

table(granulosa$coarse_annotation_archr_clean)


# trajectories ########
granulosa = addTrajectory(
  ArchRProj = granulosa, 
  name = "granulosa_pseudotime", 
  groupBy = "coarse_annotation_archr_clean",
  trajectory = c('Granulosa_sq','Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml'), 
  reducedDims = 'whole_peak_LSI',
  embedding = 'whole_peak_LSI_UMAP',
  force = TRUE
)

p = plotTrajectory(granulosa, trajectory = "granulosa_pseudotime", 
                   embedding='whole_peak_LSI_UMAP',
                   colorBy = "cellColData", name = "granulosa_pseudotime")

p[[1]]
visutils::plotVisium(umap,granulosa$granulosa_pseudotime,t='xy',xlab='',ylab='',order.points.by.z = T)

granulosa = addMotifAnnotations(ArchRProj = granulosa, motifSet = "cisbp", name = "Motif")
granulosa = addBgdPeaks(granulosa)
saveArchRProject(granulosa,outputDirectory = 'work/archr/03_granulosa_clean_v2')

# marker peaks and others ##############
granulosa = loadArchRProject('work/archr/03_granulosa_clean_v2')

# peaks = read.table('work/pycistopic/results_pycistopic_call_peaks/consensus_peaks.bed')
# peaks = GRanges(peaks$V1,ranges = IRanges(start=peaks$V2,end = peaks$V3))
# granulosa = addPeakSet(granulosa,peaks,force=TRUE)
# granulosa = addPeakMatrix(granulosa,force=TRUE)


umap = granulosa@embeddings$g_peak_UMAP$df
cex=0.7
ctcols = names(sort(table(c(granulosa$coarse_annotation_bbknn,granulosa$coarse_annotation_archr)),decreasing = TRUE))
ctcols = char2col(ctcols)[ctcols]

pdf('work/archr/03_granulosa_clean_v2/Plots/umaps.pdf',w=16,h=8)
par(mfrow=c(2,2),mar=c(1,6,1,22),xpd=NA)
visutils::plotVisium(umap,granulosa$coarse_annotation_archr,z2col=ctcols,t='xy',label.clusters = T,xlab='',ylab='',main='ArchR label transfer',show.cluster.sizes = T,cex=cex,order.points.by.z = T)
visutils::plotVisium(umap,granulosa$coarse_annotation_bbknn,z2col=ctcols,t='xy',label.clusters = T,xlab='',ylab='',main='scvi label transfer',show.cluster.sizes = T,cex=cex,order.points.by.z = T)
visutils::plotVisium(umap,granulosa$coarse_annotation_archr_clean,z2col=ctcols,t='xy',label.clusters = T,xlab='',ylab='',main='clean ArchR annotation',show.cluster.sizes = T,cex=cex,order.points.by.z = T)
visutils::plotVisium(umap,granulosa$granulosa_pseudotime,t='xy',xlab='',ylab='',order.points.by.z = T,main='pseudotime')
dev.off()

table(granulosa$Sample,granulosa$coarse_annotation_archr_clean)
table(granulosa$Donor,granulosa$coarse_annotation_archr_clean)
table(granulosa$coarse_annotation_archr_clean)


markerPeaks = getMarkerFeatures(
  ArchRProj = granulosa, 
  useMatrix = "PeakMatrix", 
  groupBy = "coarse_annotation_archr_clean",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = 'Granulosa_sq',
  useGroups = c('Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml'),
  testMethod = "wilcoxon"
)


markerPeaks_sq = getMarkerFeatures(
  ArchRProj = granulosa, 
  useMatrix = "PeakMatrix", 
  groupBy = "coarse_annotation_archr_clean",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = c('Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml'),
  useGroups = 'Granulosa_sq',
  testMethod = "wilcoxon"
)

markerPeaks = cbind(markerPeaks_sq,markerPeaks[,1:3])
dim(markerPeaks)


#plotMarkers(seMarker = markerPeaks,name='Granulosa_AMH_ml',cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
#plotMarkers(seMarker = markerPeaks,name='Granulosa_sq',cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 4",
  transpose = TRUE,
  scaleRows = TRUE,
  binaryClusterRows = FALSE,
  clusterCols = TRUE,
  nLabel = 10
)

heatmapPeaks
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", 
        width = 8, height = 6, ArchRProj = granulosa, addDOC = FALSE)


up = colSums(assay(markerPeaks,'Log2FC') >  1 & assay(markerPeaks,'FDR') < 0.05)
dw = colSums(assay(markerPeaks,'Log2FC') <  -1 & assay(markerPeaks,'FDR') < 0.05)

up4 = colSums(assay(markerPeaks,'Log2FC') >  4 & assay(markerPeaks,'FDR') < 0.05)
dw4 = colSums(assay(markerPeaks,'Log2FC') <  -4 & assay(markerPeaks,'FDR') < 0.05)


pdf('work/archr/03_granulosa_clean_v2/Plots/marker_peaks_cnt.pdf',w=8,h=5)
par(mfrow=c(1,2),mar=c(14,4,1,1),las=2)
barplotWithText(log10(up+1),up,main='Peaks up',ylab='log10',col='#EEEEEE',border=NA)
barplotWithText(log10(up4+1),up4,ylab='log10',col='#888888',den=20,add=T,border=NA)
barplotWithText(log10(dw),dw,main='Peaks down',ylab='log10',col='#EEEEEE',border=NA)
barplotWithText(log10(dw4),dw4,ylab='log10',col='#888888',den=20,add=T,border=NA)
legend(grconvertX(0,'nfc','user'),grconvertY(0,'nfc','user'),fill=c('#EEEEEE','#888888'),den=c(-1,20),legend = c('l2FC > 1','l2FC > 4'),xpd=NA,yjust = 0)
dev.off()


# _check AMH ##################
gene.descr[gene.descr$name=='AMH',]
peaks = rowData(markerPeaks)
f = peaks$seqnames=='chr19' & peaks$start>2248000 & peaks$end<2250000
peaks[f,]

assayNames(markerPeaks)
assay(markerPeaks[f,],'Log2FC')
assay(markerPeaks[f,],'Pval')
assay(markerPeaks[f,],'Mean')



p <- plotBrowserTrack(
  ArchRProj = granulosa, 
  groupBy = "coarse_annotation_archr_clean", 
  region = GRanges('chr19',IRanges(start = 2248440 -1e4,2249450+1e4)),
  #geneSymbol = c("AMH"),
  useGroups = c('Granulosa_sq','Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml'),
  features =  getMarkers(markerPeaks[f,], cutOff = "Log2FC<Inf", returnGR = TRUE)["Granulosa_sq"],
  upstream = 10000,
  downstream = 10000,
  normMethod = 'ReadsInTSS',
  tileSize=50,
)
grid::grid.newpage()
grid::grid.draw(p)

# _try my coverage plot ############
path = '/lustre/scratch127/cellgen/cellgeni/tickets/tic-3942/'
# sample_id to fragment file link
fragments_paths = read.csv(paste0(path,'actions/samples.csv'))
fragments_paths$fragment_file = paste0(fragments_paths$filedir,'/fragments.tsv.gz')

# celltype annotation
path = '/lustre/scratch127/cellgen/cellgeni/tickets/tic-3942/'
# sample_id to fragment file link
fragments_paths = read.csv(paste0(path,'actions/samples.csv'))
fragments_paths$fragment_file = paste0(fragments_paths$filedir,'/fragments.tsv.gz')

# celltype annotation
barcodes = as.data.frame(granulosa@cellColData)[,c('Sample','coarse_annotation_archr_clean')]
barcodes$barcode = splitSub(rownames(barcodes),'#',2)
colnames(barcodes) = c('sample_id','celltype','barcode')
# gene coordinates (longest protein coding isoform per gene based on 2020A 10x reference)
# gtf=readRDS(paste0(path,'gtf_longest_tr.rds')) # to plot longest protein coding isoform
gtf=readRDS(paste0(path,'gtf_all_exons_on_one_transc.rds')) # to plot all exons of the gene on one line
# cell*peak counts. loaded only for normalization normalization
celltype_totals = sapply(split(granulosa$ReadsInTSS,granulosa$coarse_annotation_archr_clean),sum)
celltype_totals[]=1

# peak to plot
peak = pasreCoors('chr19:2248440-2249450')
# celltypes to plots
cts = c('Granulosa_sq','Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml')
# margin width to include into plot
margin = 1e4
c1=getCoverage(fragments_paths,barcodes,peak,dedupl=FALSE,margin=margin)

plotCoverages(c1[cts],celltype_totals[cts]/1e6,ylab='CPM',
              region2mark = peak,
              gtf=gtf,
              xaxt='n',xlab='',ylim=NULL)

# motif enrichemnt ############
motifsUp <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = granulosa,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

motifsDo <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = granulosa,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)


sort(rownames(motifsUp))
assayNames(motifsUp)
colSums(assay(motifsUp,'mlog10Padj') >= -log10(0.0005))
colSums(assay(motifsDo,'mlog10Padj') >= -log10(0.05))

heatmapEMu <- plotEnrichHeatmap(motifsUp[,2:3], n = 10,
                                transpose = TRUE,cutOff = -log10(0.05))

heatmapEMd <- plotEnrichHeatmap(motifsDo, n = 17,
                                transpose = TRUE,cutOff = -log10(0.05))

plotPDF(heatmapEMu, name = "TF_in_up_peaks", 
        width = 8, height = 6, ArchRProj = granulosa, addDOC = FALSE)
plotPDF(heatmapEMd, name = "TF_in_down_peaks", 
        width = 8, height = 6, ArchRProj = granulosa, addDOC = FALSE)

# saveRDS(markerPeaks,'work/archr/03_granulosa_clean_v2/markerPeaks.rds')
# saveRDS(motifsUp,'work/archr/03_granulosa_clean_v2/motifsUp.rds')
# saveRDS(motifsDo,'work/archr/03_granulosa_clean_v2/motifsDo.rds')

# compare with scanpy
sc = read.csv('actions/ovary/pycistopic_granulosa_vs-sq_all-cells_all_peaks.csv',row.names = 1)
sc[1:2,]
table(sc$group)
compTests = function(sc,arch,group){
  sc = sc[sc$group==group ,]
  peaks = rowData(arch)
  rownames(arch) = paste0(peaks$seqnames,':',peaks$start,'-',peaks$end)
  arch=arch[sc$names,group]
  for(a in assayNames(arch))
    sc[[paste0('archr_',a)]] = assay(arch,a)[,1]
  sc
}

g_tr = compTests(sc,markerPeaks,'Granulosa_sq_transitioning')
table(sc=g_tr$pvals_adj<0.005,arch=g_tr$archr_FDR<0.005)
f = g_tr$pvals_adj<0.05 | g_tr$archr_FDR < 0.05
cor(g_tr$logfoldchanges,g_tr$archr_Log2FC,m='sp')
plotLine(g_tr$logfoldchanges[f],g_tr$archr_Log2FC[f],pch='.',xlab='scanpy',ylab='archr')



#  Encode TFBS #####
granulosa = addArchRAnnotations(ArchRProj = granulosa, 
                                collection = "EncodeTFBS")

enrichEncodeUp <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = granulosa,
  peakAnnotation = "EncodeTFBS",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

enrichEncodeDo <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = granulosa,
  peakAnnotation = "EncodeTFBS",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)

plotEnrichHeatmap(enrichEncodeUp, n = 7, transpose = TRUE,cutOff = -log10(0.05))
plotEnrichHeatmap(enrichEncodeDo, n = 7, transpose = TRUE,cutOff = -log10(0.05))

granulosa = addArchRAnnotations(ArchRProj = granulosa, 
                                collection = "ATAC")

enrichATACup <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = granulosa,
  peakAnnotation = "ATAC",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

enrichATACdo <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = granulosa,
  peakAnnotation = "ATAC",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)

plotEnrichHeatmap(enrichATACup, n = 7, transpose = TRUE,cutOff = -log10(0.05))
plotEnrichHeatmap(enrichATACdo, n = 7, transpose = TRUE,cutOff = -log10(0.05))

granulosa = addArchRAnnotations(ArchRProj = granulosa, collection = "Codex")

enrichCodexUp <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = granulosa,
  peakAnnotation = "Codex",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

enrichCodexDo <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = granulosa,
  peakAnnotation = "Codex",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)

plotEnrichHeatmap(enrichCodexUp, n = 7, transpose = TRUE,cutOff = -log10(0.05))
plotEnrichHeatmap(enrichCodexDo, n = 7, transpose = TRUE,cutOff = -log10(0.05))


# downsample ##########################
# marker peaks and others ##############
# granulosa = loadArchRProject('work/archr/03_granulosa_clean_v2')
# 
# downSample = function(ann,size=min(table(ann)),seed=1234){
#   set.seed(seed)
#   inx = seq_along(ann)
#   inx = unlist(lapply(split(inx,ann),sample,size=size))
#   inx
# }
# 
# inx=downSample(granulosa$coarse_annotation_archr_clean)
# 
# granulosa_sub <- subsetArchRProject(
#   ArchRProj = granulosa,
#   cells = granulosa$cellNames[inx],
#   outputDirectory = "work/archr/03_granulosa_clean_v2_sub",
#   dropCells = TRUE,
#   force = TRUE
# )
# 
# rm(granulosa)
# _check peak matrix ##################
granulosa_sub = loadArchRProject('work/archr/03_granulosa_clean_v2_sub')
amh = gsm['AMH',granulosa_sub$cellNames]
ord = c('Granulosa_sq','Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml')

pdf('work/archr/03_granulosa_clean_v2_sub/Plots/AMH_expression.pdf',w=9,h=8)
par(mfrow=c(2,2),mar=c(0,0,1,10))
plotVisium(granulosa_sub@embeddings$g_peak_UMAP$df,granulosa_sub$coarse_annotation_archr_clean,t='xy')
plotVisium(granulosa_sub@embeddings$g_peak_UMAP$df,amh,t='xy',zfun=log1p,cex=0.6,order.points.by.z = T)
par(mar=c(13,4,1,0),las=3)
boxplot(split(log1p(amh), granulosa_sub$coarse_annotation_archr_clean)[ord],ylab='log1p AMH cpm')
barplot(tapply(amh, granulosa_sub$coarse_annotation_archr_clean,mean)[ord],ylab='mean AMH cpm')
dev.off()

getAvailableMatrices(granulosa_sub)
mtx = getMatrixFromProject(granulosa_sub,'PeakMatrix')
peak = 'chr19:2248950-2249450'
peak_cnts = readRDS('work/pycistopic/results_pycistopic_call_peaks/peak_chr19_2248950-2249450_cov.rds')
i=findOverlaps(GRanges('chr19',IRanges(2248950 ,2249450)),rowRanges(mtx),type = 'any')@to
cmn = intersect(colnames(mtx),names(peak_cnts))
arch_cnts = assay(mtx,'PeakMatrix')[i,cmn]

# complete bullshit!
table(pycnt=peak_cnts[cmn],arch_cnts)

# lets add it again!
# peaks = read.table('work/pycistopic/results_pycistopic_call_peaks/consensus_peaks.bed')
# peaks = GRanges(peaks$V1,ranges = IRanges(start=peaks$V2,end = peaks$V3))
# granulosa_sub = addPeakSet(granulosa_sub,peaks,force = TRUE)
# granulosa_sub = addPeakMatrix(granulosa_sub,force=TRUE)


umap = granulosa_sub@embeddings$g_peak_UMAP$df
cex=0.7
ctcols = names(sort(table(c(granulosa_sub$coarse_annotation_bbknn,granulosa_sub$coarse_annotation_archr)),decreasing = TRUE))
ctcols = char2col(ctcols)[ctcols]

pdf('work/archr/03_granulosa_clean_v2_sub/Plots/umaps.pdf',w=16,h=8)
par(mfrow=c(2,2),mar=c(1,6,1,22),xpd=NA)
visutils::plotVisium(umap,granulosa_sub$coarse_annotation_archr,z2col=ctcols,t='xy',label.clusters = T,xlab='',ylab='',main='ArchR label transfer',show.cluster.sizes = T,cex=cex,order.points.by.z = T)
visutils::plotVisium(umap,granulosa_sub$coarse_annotation_bbknn,z2col=ctcols,t='xy',label.clusters = T,xlab='',ylab='',main='scvi label transfer',show.cluster.sizes = T,cex=cex,order.points.by.z = T)
visutils::plotVisium(umap,granulosa_sub$coarse_annotation_archr_clean,z2col=ctcols,t='xy',label.clusters = T,xlab='',ylab='',main='clean ArchR annotation',show.cluster.sizes = T,cex=cex,order.points.by.z = T)
visutils::plotVisium(umap,granulosa_sub$granulosa_pseudotime,t='xy',xlab='',ylab='',order.points.by.z = T,main='pseudotime')
dev.off()

table(granulosa_sub$Sample,granulosa$coarse_annotation_archr_clean)
table(granulosa_sub$Donor,granulosa$coarse_annotation_archr_clean)
table(granulosa_sub$coarse_annotation_archr_clean)


markerPeaks = getMarkerFeatures(
  ArchRProj = granulosa_sub, 
  useMatrix = "PeakMatrix", 
  groupBy = "coarse_annotation_archr_clean",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = 'Granulosa_sq',
  useGroups = c('Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml'),
  testMethod = "wilcoxon"
)


markerPeaks_sq = getMarkerFeatures(
  ArchRProj = granulosa_sub, 
  useMatrix = "PeakMatrix", 
  groupBy = "coarse_annotation_archr_clean",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = c('Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml'),
  useGroups = 'Granulosa_sq',
  testMethod = "wilcoxon"
)

markerPeaks = cbind(markerPeaks_sq,markerPeaks[,1:3])
dim(markerPeaks)

# what are NAs in FDR?
# seems like these are peaks with zero expression in the group
assayNames(markerPeaks)
fdr = assay(markerPeaks,'FDR')
l2fc = assay(markerPeaks,'Log2FC')
mean = assay(markerPeaks,'Mean')
table(l2fc[is.na(fdr)])
table(mean[is.na(fdr)])


#plotMarkers(seMarker = markerPeaks,name='Granulosa_AMH_ml',cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
#plotMarkers(seMarker = markerPeaks,name='Granulosa_sq',cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 4",
  transpose = TRUE,
  scaleRows = TRUE,
  binaryClusterRows = FALSE,
  clusterCols = TRUE,
  nLabel = 20
)

heatmapPeaks
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", 
        width = 8, height = 6, ArchRProj = granulosa_sub, addDOC = FALSE)


up = colSums(assay(markerPeaks,'Log2FC') >  1 & assay(markerPeaks,'FDR') < 0.05)
dw = colSums(assay(markerPeaks,'Log2FC') <  -1 & assay(markerPeaks,'FDR') < 0.05)

up4 = colSums(assay(markerPeaks,'Log2FC') >  4 & assay(markerPeaks,'FDR') < 0.05)
dw4 = colSums(assay(markerPeaks,'Log2FC') <  -4 & assay(markerPeaks,'FDR') < 0.05)


pdf('work/archr/03_granulosa_clean_v2_sub//Plots/marker_peaks_cnt.pdf',w=8,h=5)
par(mfrow=c(1,2),mar=c(14,4,1,1),las=2)
barplotWithText(log10(up+1),up,main='Peaks up',ylab='log10',col='#EEEEEE',border=NA)
barplotWithText(log10(up4+1),up4,ylab='log10',col='#888888',den=20,add=T,border=NA)
barplotWithText(log10(dw),dw,main='Peaks down',ylab='log10',col='#EEEEEE',border=NA)
barplotWithText(log10(dw4),dw4,ylab='log10',col='#888888',den=20,add=T,border=NA)
legend(grconvertX(0,'nfc','user'),grconvertY(0,'nfc','user'),fill=c('#EEEEEE','#888888'),den=c(-1,20),legend = c('l2FC > 1','l2FC > 4'),xpd=NA,yjust = 0)
dev.off()


# _check AMH ###########
gene.descr[gene.descr$name=='AMH',]
peaks = rowData(markerPeaks)
f = peaks$seqnames=='chr19' & peaks$start>2248000 & peaks$end<2250000
peaks[f,]

assayNames(markerPeaks)
assay(markerPeaks[f,],'Log2FC')
assay(markerPeaks[f,],'FDR')
assay(markerPeaks[f,],'Mean')
rownames(markerPeaks) = paste0(as.character(rowData(markerPeaks)$seqnames),':',
                               rowData(markerPeaks)$start,'-', rowData(markerPeaks)$end)

# _save markers ###########
#saveRDS(markerPeaks,'work/archr/03_granulosa_clean_v2_sub/markerPeaks.rds')
markerPeaks = readRDS('work/archr/03_granulosa_clean_v2_sub/markerPeaks.rds')
peaks_with_genes = read.csv('work/pycistopic/results_pycistopic_call_peaks/consensus_peaks_with_genes.csv',row.names = 1)

peaks = rowData(markerPeaks)
mean = assay(markerPeaks,'Mean')
fdr = assay(markerPeaks,'FDR')
l2fc = assay(markerPeaks,'Log2FC')
pv = assay(markerPeaks,'Pval')

rownames(pv) = rownames(mean) = rownames(fdr) = rownames(l2fc) = paste0(peaks$seqnames,':',peaks$start,'-',peaks$end) 
peaks_with_genes = peaks_with_genes[rownames(fdr),]

colnames(fdr) = paste0('fdr_',colnames(fdr))
colnames(l2fc) = paste0('l2fc_',colnames(l2fc))
colnames(pv) = paste0('pv_',colnames(pv))
dap_sub = cbind(peaks_with_genes,fdr,l2fc,pv)
dap_sub['chr19:2248950-2249450',]
amh=dap_sub[!is.na(dap_sub$nearest_gene_name) & dap_sub$nearest_gene_name=='AMH',]
amh[order(abs(amh$nearest_gene_dist2tss)),-5:-4]
#write.csv(dap_sub,'work/archr/03_granulosa_clean_v2_sub/markerPeaks.csv')
# 
# peaks_with_genes[fdr[,2] < 0.05 & l2fc[,2] > 4,1:4]
# l2fc[fdr[,2] < 0.05 & l2fc[,2] > 4,]

# _check couple of examples ############
# on full dataset to have higher coverage
granulosa = loadArchRProject('work/archr/03_granulosa_clean_v2')

markers_sgn = lapply(colnames(fdr),function(ct){
  f = dap_sub[,ct] < 0.05
  f[is.na(f)] = FALSE
  t = dap_sub[f,]
  t[order(abs(t[,sub('^fdr_','l2fc_',ct)]),decreasing = TRUE),]
})
names(markers_sgn) = sub('fdr_','',colnames(fdr))
markers_sgn$Granulosa_sq_transitioning[1:10,-4:-5]


toGR = function(c,mar=0){
  c = strsplit(c,'[:-]')[[1]]
  GRanges(c[1],IRanges(start = as.numeric(c[2]) - mar,as.numeric(c[3]) + mar))
}


peak = 'chr1:88907279-88907779'
for(ct in names(markers_sgn)){
  sgn = markers_sgn[[ct]]
  peaks = c(rownames(sgn[sgn[,paste0('l2fc_',ct)]>0,])[1:5],rownames(sgn[sgn[,paste0('l2fc_',ct)]<0,])[1:5])
  pdf(paste0('work/archr/03_granulosa_clean_v2_sub/Plots/examples/',ct,'.pdf'),w=10,h=5)
  for(peak in peaks){
    p <- plotBrowserTrack(
      ArchRProj = granulosa, 
      groupBy = "coarse_annotation_archr_clean", 
      region = toGR(peak,2e4),
      #geneSymbol = c("AMH"),
      useGroups = c('Granulosa_sq','Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml'),
      features =  getMarkers(markerPeaks, cutOff = "abs(Log2FC)>1 & FDR<0.05", returnGR = TRUE)["Granulosa_sq_transitioning"],
      normMethod = 'ReadsInTSS',
      tileSize=50,
    )
    
     grid::grid.newpage()
     grid::grid.draw(p)
  }
  dev.off()
}

# _cluster peaks ########
sgn = apply(fdr<0.05 & abs(l2fc)>2,1,sum)>1

mtx = getMatrixFromProject(granulosa,'PeakMatrix')
pmtx = assay(mtx,'PeakMatrix')
p = rowRanges(mtx)
rownames(pmtx) = paste0(p@seqnames,':',p@ranges@start,'-',as.integer(p@ranges@start+p@ranges@width-1))
pmtx = pmtx[rownames(fdr),]
pmtx_ct = as.matrix(visutils::calcColSums(pmtx,mtx$coarse_annotation_archr_clean))



pmtx_ct = as.matrix(visutils::calcColSums(pmtx,mtx$coarse_annotation_archr_clean))
pmtx_ct = pmtx_ct[,ord]

# table(is.na(mtx$granulosa_pseudotime),mtx$coarse_annotation_archr_clean)
# f = !is.na(mtx$granulosa_pseudotime)
# bin = number2bin(mtx$granulosa_pseudotime[f],20)
# table(bin,mtx$coarse_annotation_archr_clean[f])[,ord]
# pmtx_ct = as.matrix(visutils::calcColSums(pmtx[,f],as.character(bin)))
# pmtx_ct = pmtx_ct[,as.character(1:20)]


pmtx_ct = sweep(pmtx_ct,2,colSums(pmtx_ct),'/')*1e6

mean_sgn = log10p1(as.matrix(pmtx_ct[sgn,]))
# mean_sgn = as.matrix(mean[sgn,])

dim(mean_sgn)
cr = cor(t(mean_sgn))
hc = hclust(as.dist(1-cr))
cl = visutils::renameClustsBySize(cutree(hc,k=6))
table(cl)

par(mfrow=c(2,3))
for(c in 1:max(cl)){
  boxplot((mean_sgn[cl==c,]),main=paste0('c',c,' (',sum(cl==c),')'))
  #abline(v=c(13,16,19),lty=2,col='red')
}

dap_sub[names(cl[cl==6]),-4:-5]

# _motif enrichment #####
markerPeaks = readRDS('work/archr/03_granulosa_clean_v2_sub/markerPeaks.rds')
#granulosa_sub = addMotifAnnotations(ArchRProj = granulosa_sub, motifSet = "cisbp", annoName = "Motif",force = TRUE)
#granulosa_sub = addMotifAnnotations(ArchRProj = granulosa_sub, motifSet = "JASPAR2020", annoName =  "Motif_JASPAR2020",force = TRUE)
getAvailableMatrices(granulosa_sub)
matchesc = getMatches(ArchRProj = granulosa_sub, name = "Motif")
matchesj = getMatches(ArchRProj = granulosa_sub, name = "Motif_JASPAR2020")
rownames(matchesc) = paste0(seqnames(matchesc),':', start(matchesc),'-', end(matchesc))
matchesc = matchesc[rownames(markerPeaks),]

rownames(matchesj) = paste0(seqnames(matchesj),':', start(matchesj),'-', end(matchesj))
matchesj = matchesj[rownames(markerPeaks),]
# saveRDS(matchesj,'work/archr/03_granulosa_clean_v2_sub/matches_jaspar2020.rds')
# saveRDS(matchesc,'work/archr/03_granulosa_clean_v2_sub/matches_cisbp.rds')

matchesj['chr19:2248950-2249450',]
which(assay(matches)['chr19:2248950-2249450',])

# comp with function from ArchR
# motifsUp <- peakAnnoEnrichment(
#   seMarker = markerPeaks,
#   ArchRProj = granulosa_sub,
#   peakAnnotation = "Motif",
#   cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
# )
# 
# class(motifsUp)
# assayNames(motifsUp)
# assay(motifsUp,'BackgroundFrequency')
# looks almost identical
# myMotifsUp = do.call(cbind,lapply(ord,function(ct){
#   print(ct)
#   set = assay(markerPeaks,'FDR')[,ct] < 0.1 &  assay(markerPeaks,'Log2FC')[,ct] <= -0.5
#   myPeakAnnoEnrichment(p2m = matches,set=set,bkg=TRUE,name = ct)
# }))
# plot(assay(myMotifsUp,'mlog10p')[,2] - assay(motifsUp,'mlog10p')[,2])
# table(assay(myMotifsUp,'Enrichment')== assay(motifsUp,'Enrichment'))
# 
# motifsUp = motifsUp[,colnames(myMotifsUp)]


# _motif enchr tests  #############
# __ define sets ################
motif_enrich_tests = list()
test_sets = bkg_sets = list()


# down in Granulosa_sq_transitioning
ct ="Granulosa_sq_transitioning"
name='Granulosa_sq_transitioning_down'
test_sets[[name]] = assay(markerPeaks,'FDR')[,ct] < 0.1 &  assay(markerPeaks,'Log2FC')[,ct] <= -0.5
bkg_sets[[name]] = assay(markerPeaks,'Mean')[,'Granulosa_sq']>0.01

# up in Granulosa_AMH_early
ct='Granulosa_AMH_early'
name = 'Granulosa_AMH_early_up'
test_sets[[name]] = assay(markerPeaks,'FDR')[,ct] < 0.1 &  assay(markerPeaks,'Log2FC')[,ct] >= 0.5
bkg_sets[[name]] = rep(TRUE,nrow(markerPeaks))

# up in Granulosa_AMH_ml
ct='Granulosa_AMH_ml'
name = 'Granulosa_AMH_ml_up'
test_sets[[name]] = assay(markerPeaks,'FDR')[,ct] < 0.1 &  assay(markerPeaks,'Log2FC')[,ct] >= 0.5
bkg_sets[[name]] = rep(TRUE,nrow(markerPeaks))

# Granulosa_AMH_ml-only-down vs shared
ct='Granulosa_AMH_ml'
name = 'Granulosa_AMH_ml_down_shared_vs_only'
sgn = assay(markerPeaks,'FDR')[,ord[-1]] < 0.1 &  assay(markerPeaks,'Log2FC')[,ord[-1]] <= -0.5

test_sets[[name]] = sgn[,'Granulosa_AMH_ml'] & rowSums(sgn[,1:2])>0
bkg_sets[[name]]  = sgn[,'Granulosa_AMH_ml']

sapply(test_sets,sum)
sapply(bkg_sets,sum)

# filter peaks by being open in at least one granulosa
mtx = getMatrixFromProject(granulosa_sub,'PeakMatrix')
p = as.data.frame(rowRanges(mtx))
rownames(mtx) = paste0(p$seqnames,':',p$start,'-',p$end)
mtx = mtx[rownames(markerPeaks),]
library(Matrix)
cpm = sweep(assay(mtx,'PeakMatrix'),2,colSums(assay(mtx,'PeakMatrix')),'/')*1e4
cpm =  as.matrix(visutils::calcColSums(cpm,mtx$coarse_annotation_archr_clean,mean = TRUE)[,colnames(markerPeaks)])

means = as.matrix(visutils::calcColSums(assay(mtx,'PeakMatrix'),mtx$coarse_annotation_archr_clean))[,colnames(markerPeaks)]
means = sweep(means,2,colSums(means),'/')*1e4

# check I have same means as archr
plot(means[,1],markerPeaks@assays@data$Mean[,1],pch='.')
abline(a=0,b=1,col='red')

plot(cpm[,1],markerPeaks@assays@data$Mean[,1],pch='.')
abline(a=0,b=1,col='red')


visutils::imageWithText(cor(means,markerPeaks@assays@data$Mean))
visutils::imageWithText(cor(cpm,markerPeaks@assays@data$Mean))

sgn = apply(do.call(cbind,test_sets),1,any)
max_mean_cpm = apply(cpm,1,mean)
hist(log10(0.001+max_mean_cpm ),1000)
boxplot(log10(0.001+max_mean_cpm ) ~ sgn)
min(max_mean_cpm[sgn]) # 0.005
table(max_mean_cpm >0.005, sgn=sgn)

cmn_bkg = max_mean_cpm >0.005

# __run tests ########
tests = list(
  cisbp_ind_bkg = list(),
  cisbp_cmn_bkg = list(),
  jaspar2020_ind_bkg = list(),
  jaspar2020_cmn_bkg = list()
)
for(n in names(test_sets)){
  print(n)
  tests$cisbp_ind_bkg[[n]] = myPeakAnnoEnrichment(p2m = matchesc,set=test_sets[[n]],bkg=bkg_sets[[n]],name = n)
  tests$cisbp_cmn_bkg[[n]] = myPeakAnnoEnrichment(p2m = matchesc,set=test_sets[[n]],bkg=cmn_bkg,name = n)
  tests$jaspar2020_ind_bkg[[n]] = myPeakAnnoEnrichment(p2m = matchesj,set=test_sets[[n]],bkg=bkg_sets[[n]],name = n)
  tests$jaspar2020_cmn_bkg[[n]] = myPeakAnnoEnrichment(p2m = matchesj,set=test_sets[[n]],bkg=cmn_bkg,name = n)
}

tests = lapply(tests,function(t)do.call(cbind,t))

for(n in names(tests)){
  fdr = assay(tests[[n]],'FDR')
  pv  = assay(tests[[n]],'pv')
  or  = assay(tests[[n]],'Enrichment')
  
  colnames(fdr) = paste0('fdr_',colnames(fdr))
  colnames(pv)  = paste0('pv_',colnames(pv))
  colnames(or)  = paste0('oddsRatio_',colnames(or))
  
  write.csv(cbind(pv,or,fdr),paste0('work/archr/03_granulosa_clean_v2_sub/enriched_motifs_',n,'_v2.csv'))
}
# saveRDS(tests,'work/archr/03_granulosa_clean_v2_sub/enriched_motifs_v2.rds')

# comp different bkgs
ind = as.data.frame(myGetMarkers(tests$cisbp_ind_bkg,cutOff = "FDR<0.05 & l2or > 0.5"))
cmn = as.data.frame(myGetMarkers(tests$cisbp_cmn_bkg,cutOff = "FDR<0.05 & l2or > 0.5"))

imageWithText(compareSets(split(ind$TF_name,ind$comparison),
            split(cmn$TF_name,cmn$comparison),fun = 'j')[o,o])

compareSets(split(ind$TF_name,ind$comparison),
            split(cmn$TF_name,cmn$comparison),fun = function(x,y)length(intersect(x,y)))[o,o]

compareSets(split(ind$TF_name,ind$comparison),
            split(cmn$TF_name,cmn$comparison),fun = function(x,y)length(union(x,y)))[o,o]

# why some TFs pop up in multiple comparisons
motif_enrich_tests = tests$cisbp_cmn_bkg #readRDS('work/archr/03_granulosa_clean_v2_sub/enriched_motifs_v1.rds')
tfs = as.data.frame(myGetMarkers(motif_enrich_tests,cutOff = "FDR<0.05 & l2or > 0.5"))
table(tfs$comparison)
table(tfs$TF_name)
o=colnames(motif_enrich_tests)
tf_names = split(tfs$TF_name,tfs$comparison)
compareSets(tf_names)[o,o]
intersect(tf_names$Granulosa_sq_transitioning_down,tf_names$Granulosa_AMH_early_up)
intersect(tf_names$Granulosa_sq_transitioning_down,tf_names$Granulosa_AMH_ml_up)


ct ="Granulosa_sq_transitioning"
groups = data.frame(G_sq_open = assay(markerPeaks,'Mean')[,'Granulosa_sq']>0.01,
                    G_tr_down = assay(markerPeaks,'FDR')[,ct] < 0.1 &  assay(markerPeaks,'Log2FC')[,ct] <= -0.5)

ct='Granulosa_AMH_early'
groups$G_AMH_e_up = assay(markerPeaks,'FDR')[,ct] < 0.1 &  assay(markerPeaks,'Log2FC')[,ct] >= 0.5
ct='Granulosa_AMH_ml'
groups$G_AMH_ml_up = assay(markerPeaks,'FDR')[,ct] < 0.1 &  assay(markerPeaks,'Log2FC')[,ct] >= 0.5

groups$summary = ''
groups$summary[groups$G_AMH_e_up] = 'G_AMH_e_up'
groups$summary[groups$G_AMH_ml_up] = 'G_AMH_ml_up'
groups$summary[groups$G_AMH_e_up & groups$G_AMH_ml_up] = 'G_AMH_both_up'
groups$summary[groups$G_tr_down] = paste0(groups$summary[groups$G_tr_down],'G_tr_down')
groups$summary[groups$summary==''] = 'bkg'

table(groups$summary,groups$G_sq_open)
o=c('bkg','G_tr_down','G_AMH_e_up','G_AMH_ml_up','G_AMH_both_up')

par(mfrow=c(2,2),mar=c(10,4,1,0),las=2)
f = !groups$G_sq_open
t=table(groups$summary[f],assay(matches,'matches')[f,'SMAD5_866'])[o[o %in% groups$summary[f]],]
barplotWithText(t[,'TRUE']/rowSums(t)*100,paste0(t[,'TRUE'],'/',rowSums(t)),ylab='% of peaks with motif',main='SMAD5_866; closed in G_sq',ylim=range(0,25))

f = groups$G_sq_open
t=table(groups$summary[f],assay(matches,'matches')[f,'SMAD5_866'])[o[o %in% groups$summary[f]],]
barplotWithText(t[,'TRUE']/rowSums(t)*100,paste0(t[,'TRUE'],'/',rowSums(t)),ylab='% of peaks with motif',main='SMAD5_866; open in G_sq',ylim=range(0,25))


f = !groups$G_sq_open
t=table(groups$summary[f],assay(matches,'matches')[f,'DNMT1_301'])[o[o %in% groups$summary[f]],]
barplotWithText(t[,'TRUE']/rowSums(t)*100,paste0(t[,'TRUE'],'/',rowSums(t)),ylab='% of peaks with motif',main='DNMT1_301; closed in G_sq',ylim=range(0,25))

f = groups$G_sq_open
t=table(groups$summary[f],assay(matches,'matches')[f,'DNMT1_301'])[o[o %in% groups$summary[f]],]
barplotWithText(t[,'TRUE']/rowSums(t)*100,paste0(t[,'TRUE'],'/',rowSums(t)),ylab='% of peaks with motif',main='DNMT1_301; open in G_sq',ylim=range(0,25))

# _cPG islands ############
# downloaded from https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=3213068514_wYJAAFFyL2B1d33A6fmP2Xk8rTqr&db=hg38&c=chr7&g=cpgIslandSuper
cpg = read.csv('data/cpgIslandExt.csv.gz')
peaks = rowData(markerPeaks)
peaks = GRanges(seqnames = peaks$seqnames,IRanges(start=peaks$start,end=peaks$end))
cpg_gr = GRanges(seqnames = cpg$chrom,IRanges(start=cpg$chromStart,end=cpg$chromEnd))
ovlps = findOverlaps(peaks,cpg_gr,type='any',select='all',ignore.strand=T)
peaks_in_cpg = 1:length(peaks) %in% ovlps@from

peaks_in_cpg = matrix(peaks_in_cpg,ncol=1)
rownames(peaks_in_cpg) = 
colnames(peaks_in_cpg) = 'CpG'

cpg_tests = list(
  ind_bkg = list(),
  cmn_bkg = list()
)
for(n in names(test_sets)){
  print(n)
  cpg_tests$ind_bkg[[n]] = myPeakAnnoEnrichment(p2m = peaks_in_cpg,set=test_sets[[n]],bkg=bkg_sets[[n]],name = n)
  cpg_tests$cmn_bkg[[n]] = myPeakAnnoEnrichment(p2m = peaks_in_cpg,set=test_sets[[n]],bkg=cmn_bkg,name = n)
}

cpg_tests = lapply(cpg_tests,function(t)do.call(cbind,t))
assay(cpg_tests$ind_bkg,'pv')
assay(cpg_tests$cmn_bkg,'pv')
assay(cpg_tests$ind_bkg,'BackgroundProporition')
assay(cpg_tests$cmn_bkg,'Enrichment')


i=rbind(bkg_ind=assay(cpg_tests$ind_bkg,'BackgroundProporition'),
        bkg_cmn=assay(cpg_tests$cmn_bkg,'BackgroundProporition'),
        set=assay(cpg_tests$ind_bkg,'CompareProportion'))


par(mfrow=c(1,1),las=2,mar=c(10,4,1,1))
names = c('G_tr_down','G_AMG_early_up','G_AMG_ml_up','G_AMG_down_shared')
b=barplot(i,beside = T,names.arg = names,legend.text = c('individ bkg','common bkg','DARs'),ylim=c(0,0.35),ylab='CpG island proportion')
pv = paste0('pv=',format(assay(cpg_tests$ind_bkg,'pv'),scientific=T,digits=1))
text(b[1,],i[1,]+0.01,pv,adj=c(0,0),xpd=NA,srt=90)

pv = paste0('pv=',format(assay(cpg_tests$cmn_bkg,'pv'),scientific=T,digits=1))
text(b[2,],i[2,]+0.01,pv,adj=c(0,0),xpd=NA,srt=90)


barplot(c,beside = T,names.arg = names,legend.text = T,main='Common backgrounds (open in at least one G)',ylim=c(0,0.35))

text(colSums(b)/2,apply(c,2,max)+0.01,pv,adj=c(0.5,0),xpd=NA)


lapply(test_sets,function(x)table(x,peaks_in_cpg[,1]))
lapply(bkg_sets,function(x)table(x,peaks_in_cpg[,1]))

# _cromvar ##############
granulosa_sub@peakAnnotation

granulosa_sub = addMotifAnnotations(ArchRProj = granulosa_sub, motifSet = "cisbp", name = "Motif",force = TRUE)
granulosa_sub = addBgdPeaks(granulosa_sub,force = TRUE)

granulosa_sub = addDeviationsMatrix(
  ArchRProj = granulosa_sub, 
  peakAnnotation = "Motif",
  force = TRUE
)



dev_plot = getVarDeviations(granulosa_sub,plot=TRUE)
dev_plot
plotPDF(dev_plot, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = granulosa_sub, addDOC = FALSE)

devdf = getVarDeviations(granulosa_sub,plot=FALSE)
devdf[1:10,]



p <- plotGroups(ArchRProj = granulosa_sub, 
                groupBy = "coarse_annotation_archr_clean", 
                colorBy = "MotifMatrix", 
                name = c('z:JUNB_139'),
                imputeWeights = NULL#,getImputeWeights(granulosa_sub)
)
p

p <- plotEmbedding(
  ArchRProj = granulosa_sub, 
  colorBy = "MotifMatrix", 
  name = c('z:FOXL1_369'), 
  embedding = "g_peak_UMAP",
  imputeWeights = NULL#getImputeWeights(projHeme5)
)
p


getAvailableMatrices(granulosa_sub)

markerMotifs = getMarkerFeatures(
  ArchRProj = granulosa_sub, 
  useMatrix = "MotifMatrix", 
  groupBy = "coarse_annotation_archr_clean",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = 'Granulosa_sq',
  useGroups = c('Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml'),
  testMethod = "wilcoxon"
)

assayNames(markerMotifs)
apply(assay(markerMotifs,'FDR')<0.05 & abs(assay(markerMotifs,'MeanDiff')) > 0.1,2,sum)

mm = getMarkers(markerMotifs,cutOff='FDR<0.05 & MeanDiff >0')
mm$Granulosa_sq_transitioning
mm$Granulosa_AMH_ml[order(mm$Granulosa_AMH_ml$MeanDiff,decreasing = T),]

mmz = getMatrixFromProject(granulosa_sub,'MotifMatrix')
gs = getMatrixFromProject(granulosa_sub,'GeneScoreMatrix')
assayNames(mmz)
assayNames(gs)
rownames(gs) = rowData(gs)$name
boxplot(split(assay(mmz,'z')['FOS_137',],mmz$coarse_annotation_archr_clean)[ord])
boxplot(split(log10p1(assay(gs,'GeneScoreMatrix')['ZNF32',]),mmz$coarse_annotation_archr_clean)[ord])
sapply(split(log10p1(assay(gs,'GeneScoreMatrix')['ZNF32',]),mmz$coarse_annotation_archr_clean)[ord],mean)

ge = gex$X['ZNF32',]/colSums(gex$X)*1e4
boxplot(split(log10p1(ge),gex$obs$coarse_annotation)[ord])
tapply(ge>0,gex$obs$coarse_annotation,mean)[ord]
table(gex$obs$coarse_annotation,ge>0)[ord,]

# _subset peaks to expressed ones ##############
# as it should affect background
granulosa_sub = loadArchRProject('work/archr/03_granulosa_clean_v2_sub')
granulosa_sub = saveArchRProject(granulosa_sub,'work/archr/03_granulosa_clean_v2_sub_peaks',load = TRUE)

mtx = getMatrixFromProject(granulosa_sub,'PeakMatrix')
ncells = rowSums(assay(mtx,'PeakMatrix')>0)
hist(ncells,100)
table(ncells>19)
# subset peaks to expressed
peaks = rowRanges(mtx)[ncells>19,]
granulosa_sub = addPeakSet(granulosa_sub,peaks,force = TRUE)
granulosa_sub = addPeakMatrix(granulosa_sub,force=TRUE)

granulosa_sub = addMotifAnnotations(ArchRProj = granulosa_sub, motifSet = "cisbp", name = "Motif",force = TRUE)
granulosa_sub = addBgdPeaks(granulosa_sub,force = TRUE)

granulosa_sub = addDeviationsMatrix(
  ArchRProj = granulosa_sub, 
  peakAnnotation = "Motif",
  force = TRUE
)

dev_plot = getVarDeviations(granulosa_sub,plot=TRUE)
dev_plot
plotPDF(dev_plot, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = granulosa_sub, addDOC = FALSE)

devdf = getVarDeviations(granulosa_sub,plot=FALSE)
devdf[1:10,]
devdf[devdf$name=='FOXL1_369',]



p <- plotGroups(ArchRProj = granulosa_sub, 
                groupBy = "coarse_annotation_archr_clean", 
                colorBy = "MotifMatrix", 
                name = c('z:ZNF32_237'),
                groupOrder=ord,
                imputeWeights = NULL#,getImputeWeights(granulosa_sub)
)
p

p <- plotEmbedding(
  ArchRProj = granulosa_sub, 
  colorBy = "MotifMatrix", 
  name = c('z:ZNF32_237'), 
  embedding = "g_peak_UMAP",
  imputeWeights = NULL#getImputeWeights(projHeme5)
)



getAvailableMatrices(granulosa_sub)

markerMotifs = getMarkerFeatures(
  ArchRProj = granulosa_sub, 
  useMatrix = "MotifMatrix", 
  groupBy = "coarse_annotation_archr_clean",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = 'Granulosa_sq',
  useGroups = c('Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml'),
  testMethod = "wilcoxon"
)

assayNames(markerMotifs)
apply(assay(markerMotifs,'FDR')<0.05 & abs(assay(markerMotifs,'MeanDiff')) > 0.1,2,sum)

mm = getMarkers(markerMotifs,cutOff='FDR<0.05 & MeanDiff >0')
mm$Granulosa_sq_transitioning
mm$Granulosa_AMH_ml[order(mm$Granulosa_AMH_ml$MeanDiff,decreasing = T),]
mm$Granulosa_AMH_ml[mm$Granulosa_AMH_ml$name=='FOXL1_369',]



mmz = getMatrixFromProject(granulosa_sub,'MotifMatrix')
gs = getMatrixFromProject(granulosa_sub,'GeneScoreMatrix')
assayNames(mmz)
assayNames(gs)
rownames(gs) = rowData(gs)$name
boxplot(split(assay(mmz,'z')['ZNF32_237',],mmz$coarse_annotation_archr_clean)[ord])
boxplot(split(log10p1(assay(gs,'GeneScoreMatrix')['ZNF32',]),mmz$coarse_annotation_archr_clean)[ord])
sapply(split(log10p1(assay(gs,'GeneScoreMatrix')['ZNF32',]),mmz$coarse_annotation_archr_clean)[ord],mean)

ge = gex$X['ZNF32',]/colSums(gex$X)*1e4
boxplot(split(log10p1(ge),gex$obs$coarse_annotation)[ord])
tapply(ge>0,gex$obs$coarse_annotation,mean)[ord]
table(gex$obs$coarse_annotation,ge>0)[ord,]




# motif enrichemnt ############
motifsUp <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = granulosa,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

motifsDo <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = granulosa,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)


sort(rownames(motifsUp))
assayNames(motifsUp)
colSums(assay(motifsUp,'mlog10Padj') >= -log10(0.0005))
colSums(assay(motifsDo,'mlog10Padj') >= -log10(0.05))

heatmapEMu <- plotEnrichHeatmap(motifsUp[,2:3], n = 10,
                                transpose = TRUE,cutOff = -log10(0.05))

heatmapEMd <- plotEnrichHeatmap(motifsDo, n = 17,
                                transpose = TRUE,cutOff = -log10(0.05))

plotPDF(heatmapEMu, name = "TF_in_up_peaks", 
        width = 8, height = 6, ArchRProj = granulosa, addDOC = FALSE)
plotPDF(heatmapEMd, name = "TF_in_down_peaks", 
        width = 8, height = 6, ArchRProj = granulosa, addDOC = FALSE)

# saveRDS(markerPeaks,'work/archr/03_granulosa_clean_v2/markerPeaks.rds')
# saveRDS(motifsUp,'work/archr/03_granulosa_clean_v2/motifsUp.rds')
# saveRDS(motifsDo,'work/archr/03_granulosa_clean_v2/motifsDo.rds')
write.csv(assay(motifsUp,'Enrichment'),'work/archr/03_granulosa_clean_v2/motifsUp_Enrichment.csv')
write.csv(assay(motifsUp,'mlog10Padj'),'work/archr/03_granulosa_clean_v2/motifsUp_mlog10Padj.csv')

# compare with scanpy
sc = read.csv('actions/ovary/pycistopic_granulosa_vs-sq_all-cells_all_peaks.csv',row.names = 1)
sc[1:2,]
table(sc$group)
compTests = function(sc,arch,group){
  sc = sc[sc$group==group ,]
  peaks = rowData(arch)
  rownames(arch) = paste0(peaks$seqnames,':',peaks$start,'-',peaks$end)
  arch=arch[sc$names,group]
  for(a in assayNames(arch))
    sc[[paste0('archr_',a)]] = assay(arch,a)[,1]
  sc
}

g_tr = compTests(sc,markerPeaks,'Granulosa_sq_transitioning')
table(sc=g_tr$pvals_adj<0.005,arch=g_tr$archr_FDR<0.005)
f = g_tr$pvals_adj<0.05 | g_tr$archr_FDR < 0.05
cor(g_tr$logfoldchanges,g_tr$archr_Log2FC,m='sp')
plotLine(g_tr$logfoldchanges[f],g_tr$archr_Log2FC[f],pch='.',xlab='scanpy',ylab='archr')



