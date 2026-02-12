library(ArchR)
library(Seurat)
library(visutils)
source('actions/ovary/bin/archr_utils.R')

addArchRGenome("hg38")
addArchRThreads(threads = 10)

granulosa_celltypes = c('Granulosa_sq','Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml')
gsm = readRDS('work/archr_with_fetal_atac/02_annotate/GeneScoreMatrix.rds')
gsm = unwrapAssay(gsm)
#pm = readRDS('work/archr2/03_per_class_annotation/Granulosa/peaks2cell.rds')
#prj = loadArchRProject('work/archr2/03_per_class_annotation/Granulosa/')
prj = loadArchRProject('work/archr2/04_per_class_clean/Granulosa')
#pm = getMatrixFromProject(prj,'PeakMatrix')
gene.descr = readRDS('/nfs/cellgeni/pasham/code/nf-scsajr/ref/human_2020A/functional_annotation/gene.descr.rds')

granulosa_sub = loadArchRProject(path="work/archr2/05_per_class_clean_sub/Granulosa")
#saveArchRProject(granulosa_sub,outputDirectory = "work/archr2/05_per_class_clean_sub/Granulosa")

pdf('work/archr2/03_per_class_annotation/Granulosa/figures/umaps.pdf',w=5*11,h=4*4)
ctcols = char2col(c(prj$class.coarse_annotation_archr,prj$class.coarse_annotation_archr_h))
par(mfrow=c(4,11),mar=c(0,0,1,16))
for(e in names(prj@embeddings)){
  umap = prj@embeddings[[e]]$df
  cl = sub('_UMAP','_cl_2',e)
  if(!(cl %in% colnames(prj@cellColData)))
    cl = 'TileMatrix_LSI_harmony_cl_2'
  
  plotVisium(umap,prj$Donor,t='xy',main=e)
  plotVisium(umap,prj$nFrags,main='nFrags',t='xy')
  plotVisium(umap,prj$nFrags,main='nFrags',t='xy',zfun = log1p)
  plotVisium(umap,prj@cellColData[,cl],main=cl,t='xy',zfun = log1p,label.clusters = T)
  plotVisium(umap,prj$DoubletScore,main='DoubletScore',t='xy')
  plotVisium(umap,prj$is_doublet,main='is_doublet',t='xy')
  plotVisium(umap,prj$TSSEnrichment,main='TSSEnrichment',t='xy')
  for(a in c('class.coarse_annotation_archr','coarse_annotation_archr','class.coarse_annotation_archr_h'))
    plotVisium(umap,prj@cellColData[,a],main=a,t='xy',label.clusters = T,z2col = ctcols)
  
  plotVisium(umap,amh,t='xy',main='AMH',zfun = log1p,cex=0.8,order.points.by.z = T)
}
dev.off()


# clean by QC#############
table(prj$TileMatrix_LSI_harmony_cl_2,prj$is_doublet)
table(prj$TileMatrix_LSI_harmony_cl_2,prj$nFrags>1e4)

cl_stat = paste0('C',1:length(unique(prj$TileMatrix_LSI_harmony_cl_2)))

cl_stat = data.frame(
  doublets = tapply(prj$is_doublet,prj$TileMatrix_LSI_harmony_cl_2,sum)[cl_stat],
  good_cov = tapply(prj$nFrags>1e4,prj$TileMatrix_LSI_harmony_cl_2,sum)[cl_stat],
  N = as.numeric(table(prj$TileMatrix_LSI_harmony_cl_2)[cl_stat])
)

plot(cl_stat$doublets/cl_stat$N,cl_stat$good_cov/cl_stat$N,t='n')
text(cl_stat$doublets/cl_stat$N,cl_stat$good_cov/cl_stat$N,rownames(cl_stat))


# remove clusters with more than 50% of doublets or low coverage
cl2keep = rownames(cl_stat)[cl_stat$doublets/cl_stat$N < 0.5 & cl_stat$good_cov/cl_stat$N > 0.5]
setdiff(rownames(cl_stat),cl2keep)
umap = prj@embeddings$TileMatrix_LSI_harmony_UMAP$df

pdf('work/archr2/03_per_class_annotation/Granulosa/figures/qc_cleaning_umaps.pdf',w=4*3,h=3.5*3)
par(mfrow=c(3,3),mar=c(0,0,1,7),cex=0.6)
plotVisium(umap,prj$Donor,t='xy',main=e)
plotVisium(umap,prj$nFrags,main='nFrags',t='xy')
plotVisium(umap,prj$TileMatrix_LSI_harmony_cl_2,main='TileMatrix_LSI_harmony_cl_2',t='xy',zfun = log1p,label.clusters = T)
plotVisium(umap,prj$DoubletScore,main='DoubletScore',t='xy')
plotVisium(umap,prj$is_doublet,main='is_doublet',t='xy')
plotVisium(umap,prj$TSSEnrichment,main='TSSEnrichment',t='xy')
plotVisium(umap,prj$class.coarse_annotation_archr_h,main='class.coarse_annotation_archr_h',t='xy')
plotVisium(umap,prj$TileMatrix_LSI_harmony_cl_2 %in% cl2keep,t='xy',legend.args = list(title='to keep'))
dev.off()

prj = ArchR::subsetArchRProject(prj,
                                cells = prj$cellNames[prj$TileMatrix_LSI_harmony_cl_2 %in% cl2keep],
                                outputDirectory = 'work/archr2/04_per_class_clean/Granulosa',
                                dropCells = TRUE,force = TRUE)

# check peak mtx => all good(and any)
# pm_ = readRDS('work/archr2/03_per_class_annotation/Granulosa/peaks2cell.rds')
# pm = getMatrixFromProject(prj,'PeakMatrix')
# pm = unwrapAssay(pm)
# pm_ = unwrapAssay(pm_)
# compareMtxs(pm,pm_)

# clean old clusters
cl = colnames(prj@cellColData)
cl = cl[grep('Matrix',cl)]
prj@cellColData[,cl] = NULL

prj@reducedDims = SimpleList()
prj@embeddings = SimpleList()

prj = run_rdim_cl(prj,batch='Donor',mtx_names=c('TileMatrix','PeakMatrix'))


# update gene score matrix to new ann
# geneAnnotation = myCreateGeneAnnotation('hg38','/software/cellgen/cellgeni/refdata_10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf')
# 
# prj = addGeneScoreMatrix(prj,genes = geneAnnotation$genes,force=TRUE)
# 
gsm = getMatrixFromProject(prj,'GeneScoreMatrix')
dim(gsm)

saveArchRProject(prj,outputDirectory = 'work/archr2/04_per_class_clean/Granulosa')

# annotate one more time ###########
rna = schard::h5ad2seurat('data/ref_rawcounts.h5ad')
rna$class = visutils::splitSub(rna$coarse_annotation,'_',1)
rna$class[rna$class %in% c('Theca','Stroma','FibC3')] = 'Mesenchymal'


rna_sub = rna[,rna$class=='Granulosa']
# for ngsa
rna_sub = rna_sub[,rna_sub$coarse_annotation!='Granulosa_sq_atr']
# for g4
rna_sub = rna_sub[,rna_sub$coarse_annotation %in% granulosa_celltypes]
table(rna_sub$coarse_annotation)

prj <- addGeneIntegrationMatrix(
  ArchRProj = prj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "TileMatrix_LSI_harmony",
  seRNA = rna_sub,
  addToArrow = FALSE,
  force = TRUE,
  groupRNA = 'coarse_annotation',
  nameCell = "class-coarse_annotation_g4_archr_h_cell",
  nameGroup = "class-coarse_annotation_g4_archr_h",
  nameScore = "class-coarse_annotation_g4_archr_h_score"
)


prj <- addGeneIntegrationMatrix(
  ArchRProj = prj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "TileMatrix_LSI",
  seRNA = rna_sub,
  addToArrow = FALSE,
  force = TRUE,
  groupRNA = 'coarse_annotation',
  nameCell = "class-coarse_annotation_g4_archr_cell",
  nameGroup = "class-coarse_annotation_g4_archr",
  nameScore = "class-coarse_annotation_g4_archr_score"
)

# check old ann 
p = loadArchRProject('work/archr/03_granulosa_clean_v2')
plotVisium(p@embeddings$g_peak_UMAP$df,
           p$coarse_annotation_archr_clean,main=a,t='xy',label.clusters = T,z2col = ctcols)

prj@cellColData[,'coarse_annotation_archr_clean_old'] = 'NA'
cmn = intersect(prj$cellNames,p$cellNames)  
prj@cellColData[cmn,'coarse_annotation_archr_clean_old'] = p@cellColData[cmn,'coarse_annotation_archr_clean']



amh=agm['AMH',prj$cellNames]

pdf('work/archr2/04_per_class_clean/Granulosa/figures/umaps.pdf',w=5*15,h=4*4)
ctcols = char2col(c(prj$class.coarse_annotation_archr,prj$class.coarse_annotation_archr_h))
par(mfrow=c(4,15),mar=c(0,0,1,16))
for(e in names(prj@embeddings)){
  umap = prj@embeddings[[e]]$df
  cl = sub('_UMAP','_cl_10',e)
  if(!(cl %in% colnames(prj@cellColData)))
    cl = 'TileMatrix_LSI_harmony_cl_2'
  
  plotVisium(umap,prj$Donor,t='xy',main=e)
  plotVisium(umap,prj$nFrags,main='nFrags',t='xy')
  plotVisium(umap,prj$nFrags,main='nFrags',t='xy',zfun = log1p)
  plotVisium(umap,prj@cellColData[,cl],main=cl,t='xy',zfun = log1p,label.clusters = T)
  plotVisium(umap,prj$DoubletScore,main='DoubletScore',t='xy')
  plotVisium(umap,prj$is_doublet,main='is_doublet',t='xy')
  plotVisium(umap,prj$TSSEnrichment,main='TSSEnrichment',t='xy')
  for(a in c('coarse_annotation_archr','coarse_annotation_archr_h','class.coarse_annotation_archr','class.coarse_annotation_archr_h',
             "class-coarse_annotation_g4_archr_h","class-coarse_annotation_g4_archr",
             "coarse_annotation_archr_clean_old"))
    plotVisium(umap,prj@cellColData[,a],main=a,t='xy',label.clusters = T,z2col = ctcols)
  
  plotVisium(umap,amh,t='xy',main='AMH',zfun = log1p,cex=0.8,order.points.by.z = T)
}
dev.off()


# C14 is AUrv18 and maybe Granulosa_AMH_antral_HedgLow (according to some of label transfers)
table(prj$class.coarse_annotation_archr,prj$TileMatrix_LSI_harmony_cl_10)
table(prj$Donor,prj$TileMatrix_LSI_harmony_cl_10)

# rename C14 to Granulosa_AMH_antral_HedgLow
# add old transitioning
prj@cellColData[,'class.coarse_annotation_g4_archr_h_'] = prj@cellColData[,'class.coarse_annotation_g4_archr_h']
prj@cellColData[prj$TileMatrix_LSI_harmony_cl_10=='C14','class.coarse_annotation_g4_archr_h_'] = 'Granulosa_AMH_antral_HedgLow'
prj@cellColData[prj$coarse_annotation_archr_clean_old=='Granulosa_sq_transitioning','class.coarse_annotation_g4_archr_h_'] = 'Granulosa_sq_transitioning'

prj@cellColData[,'coarse_annotation_final'] = getMajorCelltype(prj$TileMatrix_LSI_harmony_cl_1,prj@cellColData[,'class.coarse_annotation_g4_archr_h_'])



pdf('work/archr2/04_per_class_clean/Granulosa/figures/annotation_umaps.pdf',w=6*2+1,h=3.5*2)
e='TileMatrix_LSI_harmony_UMAP'
umap = prj@embeddings[[e]]$df
par(mfrow=c(2,2),mar=c(0,0,1,18),cex=0.6,oma=c(0,6,0,0),xpd=NA)
plotVisium(umap,prj$Donor,main='label transfer',t='xy',label.clusters = T,ylab='')
plotVisium(umap,prj$class.coarse_annotation_g4_archr_h,main='label transfer',t='xy',label.clusters = T,z2col = ctcols[granulosa_celltypes],show.cluster.sizes = T,ylab='')
plotVisium(umap,prj$coarse_annotation_final,main='coarse_annotation_final',t='xy',label.clusters = T,
           z2col = ctcols[c(granulosa_celltypes,'Granulosa_AMH_antral_HedgLow')],show.cluster.sizes = T,ylab='')
plotVisium(umap,prj$coarse_annotation_archr_clean_old,main='coarse_annotation_archr_clean_old',t='xy',label.clusters = T,
           z2col = c(ctcols[granulosa_celltypes],'NA'='black'),show.cluster.sizes = T,ylab='')
dev.off()

t =table(old=prj$coarse_annotation_archr_clean_old,prj$coarse_annotation_final)[c(granulosa_celltypes,'NA'),c(granulosa_celltypes,'Granulosa_AMH_antral_HedgLow')]
par(mar=c(12,14,1,1),las=3,mgp=c(11,0.3,0),tcl=-0.2)
imageWithText(log(1+t),t,ylab='new',xlab='old')

plotVisium(umap,prj$Donor=='AUrv18',main='label transfer',t='xy',label.clusters = T)

# marker genes #######
markersGSgra <- getMarkerFeatures(
  ArchRProj = prj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "coarse_annotation_final",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

#saveRDS(markersGSgra,'work/archr2/04_per_class_clean/Granulosa/markersGSgra.rds')

rownames(markersGSgra) = make.unique(rowData(markersGSgra)$name)
top5 = myGetMarkers(markersGSgra, cutOff = "FDR <= 0.05 & Log2FC >= 0.6",n=10)

table(top5$comparison)

top5l = split(top5$name,top5$comparison)
top5l$Granulosa_AMH_early = c(top5l$Granulosa_AMH_early,'AMH')

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGSgra[unlist(top5l[granulosa_celltypes]),c(granulosa_celltypes,'Granulosa_AMH_antral_HedgLow')], 
  cutOff = "FDR <= 10", 
  transpose = TRUE,clusterCols=FALSE,nLabel = 1000
)
heatmapGS


# _gene from Luz ###########
markers_granulosa = list(
  'Granulosa_sq' = c("NOTCH3","RDH10", 'COL4A4', 'HEYL','HES4'),
  'Granulosa_sq_transitioning'= c('AMH', 'DHH', 'FSHR', 'FST', 'NECTIN1', 'TNNI3', 'GJA1', "PRAME", "NEDD9"),
  'Granulosa_AMH_early'= c('MYL7', 'GSTA1', 'LRP8',  'COL2A1', 'NKAIN2'),
  'Granulosa_AMH_ml'= c('HSD11B2',  "CYP19A1",  'IHH',  'INHBB',  'MMP11', 'VCAN', 'CHGB', 'PLP1', 'IL17RA', 'TLL2'),
  'Granulosa antral'= c( 'CRHBP', 'GDF7',  'HES6',  'ZNF774', 'CXCR4', 'HSD3B2'),
  'Granulosa cummulus'= c('WFDC6',  'EPPIN',  'TTYH1')
)

markers_granulosao = names(markers_granulosa)
markers_granulosa = unlist(markers_granulosa)
markers_granulosa = data.frame(celltype=names(markers_granulosa),marker=markers_granulosa)
markers_granulosa$celltype = sub("\\d+$",'',markers_granulosa$celltype)


setdiff(markers_granulosa$marker,rownames(markersGSgra))
cols=list(celltype=char2col(markers_granulosa$celltype)[markers_granulosao])


heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGSgra[markers_granulosa$marker,c(granulosa_celltypes,'Granulosa_AMH_antral_HedgLow')], 
  cutOff = "FDR <= 10", 
  transpose = TRUE,clusterCols=FALSE,nLabel = 1000
)
heatmapGS = HeatmapAnnotation(celltype=markers_granulosa$celltype,col=cols,which='column')  %v% heatmapGS  
heatmapGS

# check donors ####

ad2c = table(prj$Donor,prj$coarse_annotation_final)[,granulosa_celltypes]
par(mar=c(4,15,1,1),bty='n')
#imageWithText(sweep(ad2c,1,rowSums(ad2c),'/'),ad2c)
imageWithText(log1p(sweep(ad2c,2,colSums(ad2c),'/')*10),ad2c)

rd2c = table(rna_sub$Donor,rna_sub$coarse_annotation)[,granulosa_celltypes]

cmnd = intersect(rownames(ad2c),rownames(rd2c))
ad2c[setdiff(rownames(ad2c),cmnd),]
setdiff(rownames(rd2c),cmnd)

ad2c = ad2c[cmnd,]
rd2c = rd2c[cmnd,]
imageWithText(log1p(sweep(ad2c,2,colSums(ad2c),'/')*10),ad2c)
imageWithText(log1p(sweep(rd2c,2,colSums(rd2c),'/')*10),rd2c)

at = rowSums(ad2c)
rt = rowSums(rd2c)
par(mfrow=c(2,2),mar=c(4,4,1,1),bty='n')
for(ct in granulosa_celltypes){
  plotLine(ad2c[,ct]/at,rd2c[,ct]/rt,main=ct,pch=16,xlab='ATAC cell fraction',ylab='RNA cell fraction',plot.ci = T)
}

saveArchRProject(prj,outputDirectory = 'work/archr2/04_per_class_clean/Granulosa')

# downsample ##########################
# granulosa = loadArchRProject('work/archr2/04_per_class_clean/Granulosa')
# 
# downSample = function(ann,size=min(table(ann)),seed=1234){
#   set.seed(seed)
#   inx = seq_along(ann)
#   inx = unlist(lapply(split(inx,ann),sample,size=size))
#   inx
# }
# 
# cells = granulosa@cellColData
# cells = cells[cells$coarse_annotation_final!='Granulosa_AMH_antral_HedgLow',]
# 
# inx=downSample(cells$coarse_annotation_final)
# table(cells$coarse_annotation_final[inx])
# 
# granulosa_sub = subsetArchRProject(
#   ArchRProj = granulosa,
#   cells = rownames(cells)[inx],
#   outputDirectory = "work/archr2/05_per_class_clean_sub/Granulosa",
#   dropCells = TRUE,
#   force = TRUE
# )

granulosa_sub = loadArchRProject(path="work/archr2/05_per_class_clean_sub/Granulosa")
table(granulosa_sub$coarse_annotation_final)

# plot AMH expression ######
amh = gsm['AMH',granulosa_sub$cellNames]
ord = c('Granulosa_sq','Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml')

pdf('work/archr2/05_per_class_clean_sub/Granulosa/figures/AMH_expression.pdf',w=11,h=8)
par(mfrow=c(2,2),mar=c(0,0,1,10),oma=c(1,1,1,1))
plotVisium(granulosa_sub@embeddings$TileMatrix_LSI_harmony_UMAP$df,granulosa_sub$coarse_annotation_final,t='xy',show.cluster.sizes = T)
plotVisium(granulosa_sub@embeddings$TileMatrix_LSI_harmony_UMAP$df,amh,t='xy',zfun=log1p,cex=0.6,order.points.by.z = T,legend.args = list(title='AMH cpm'))
par(mar=c(13,4,1,0),las=3)
boxplot(split(log1p(amh), granulosa_sub$coarse_annotation_final)[ord],ylab='log1p AMH cpm')
barplot(tapply(amh, granulosa_sub$coarse_annotation_final,mean)[ord],ylab='mean AMH cpm')
dev.off()


# _check peak matrix ##################
getAvailableMatrices(granulosa_sub)
mtx = getMatrixFromProject(granulosa_sub,'PeakMatrix')
mtx_ = unwrapAssay(mtx)

gene.descr[gene.descr$name=='AMH',]
peaks = as.data.frame(rowRanges(mtx))
f = peaks$seqnames=='chr19' & peaks$start>2248000 & peaks$end<2250000
peaks[f,]
table(mtx_[f,])

# chr19:2248841-2249341 vs
# chr19:2248950-2249450
peak_cnts = readRDS('work/pycistopic/results_pycistopic_call_peaks/peak_chr19_2248950-2249450_cov.rds')
cmn = intersect(colnames(mtx),names(peak_cnts))
arch_cnts = assay(mtx,'PeakMatrix')[i,cmn]

# looks ok
table(pycnt=peak_cnts[cmn],archr_new = mtx_[f,cmn])


# _marker peaks #############
markerPeaks = getMarkerFeatures(
  ArchRProj = granulosa_sub, 
  useMatrix = "PeakMatrix", 
  groupBy = "coarse_annotation_final",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = 'Granulosa_sq',
  useGroups = c('Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml'),
  testMethod = "wilcoxon"
)


markerPeaks_sq = getMarkerFeatures(
  ArchRProj = granulosa_sub, 
  useMatrix = "PeakMatrix", 
  groupBy = "coarse_annotation_final",
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


# heatmapPeaks <- plotMarkerHeatmap(
#   seMarker = markerPeaks, 
#   cutOff = "FDR <= 0.05 & Log2FC >= 4",
#   transpose = TRUE,
#   scaleRows = TRUE,
#   binaryClusterRows = FALSE,
#   clusterCols = TRUE,
#   nLabel = 20
# )
# 
# heatmapPeaks
# plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", 
#         width = 8, height = 6, ArchRProj = granulosa_sub, addDOC = FALSE)


up = colSums(assay(markerPeaks,'Log2FC') >  1 & assay(markerPeaks,'FDR') < 0.05)
dw = colSums(assay(markerPeaks,'Log2FC') <  -1 & assay(markerPeaks,'FDR') < 0.05)

up4 = colSums(assay(markerPeaks,'Log2FC') >  4 & assay(markerPeaks,'FDR') < 0.05)
dw4 = colSums(assay(markerPeaks,'Log2FC') <  -4 & assay(markerPeaks,'FDR') < 0.05)


pdf('work/archr2/05_per_class_clean_sub/Granulosa/figures/marker_peaks_cnt.pdf',w=8,h=5)
par(mfrow=c(1,2),mar=c(14,4,1,1),las=2)
barplotWithText(log10(up+1),up,main='Peaks up',ylab='log10',col='#EEEEEE',border=NA)
barplotWithText(log10(up4+1),up4,ylab='log10',col='#888888',den=20,add=T,border=NA)
barplotWithText(log10(dw),dw,main='Peaks down',ylab='log10',col='#EEEEEE',border=NA)
barplotWithText(log10(dw4),dw4,ylab='log10',col='#888888',den=20,add=T,border=NA)
legend(grconvertX(0,'nfc','user'),grconvertY(0,'nfc','user'),fill=c('#EEEEEE','#888888'),den=c(-1,20),legend = c('l2FC > 1','l2FC > 4'),xpd=NA,yjust = 0)
dev.off()


p = rowData(markerPeaks)
rownames(markerPeaks) = paste0(p$seqnames,':',p$start,'-',p$end)

# _save markers ###########
# saveRDS(markerPeaks,paste0(granulosa_sub@projectMetadata$outputDirectory,'/markerPeaks.rds'))
markerPeaks = readRDS(paste0(granulosa_sub@projectMetadata$outputDirectory,'/markerPeaks.rds'))
peaks_with_genes_ = read.csv('work/pycistopic/results_pycistopic_call_peaks/consensus_peaks_with_genes.csv',row.names = 1)

peaks = rowData(markerPeaks)
mean = assay(markerPeaks,'Mean')
fdr = assay(markerPeaks,'FDR')
l2fc = assay(markerPeaks,'Log2FC')
pv = assay(markerPeaks,'Pval')

rownames(pv) = rownames(mean) = rownames(fdr) = rownames(l2fc) = paste0(peaks$seqnames,':',peaks$start,'-',peaks$end) 
peaks_with_genes = granulosa_sub@peakSet
names(peaks_with_genes) = NULL
peaks_with_genes = as.data.frame(peaks_with_genes)
rownames(peaks_with_genes) = paste0(peaks_with_genes$seqnames,':',peaks_with_genes$start,'-',peaks_with_genes$end)

peaks_with_genes = peaks_with_genes[rownames(fdr),]


colnames(fdr) = paste0('fdr_',colnames(fdr))
colnames(l2fc) = paste0('l2fc_',colnames(l2fc))
colnames(pv) = paste0('pv_',colnames(pv))
dap_sub = cbind(peaks_with_genes,fdr,l2fc,pv)
amh=dap_sub[!is.na(dap_sub$nearestGene) & dap_sub$nearestGene=='AMH',]
amh[order(abs(amh$distToTSS)),-5:-4]
dap_sub['chr19:2248841-2249341',]
# write.csv(dap_sub,paste0(granulosa_sub@projectMetadata$outputDirectory,'/markerPeaks.csv'))


# __two flavours of Transitioning  #############
granulosa_sub$coarse_annotation_final_tr = granulosa_sub$coarse_annotation_final
granulosa_sub$coarse_annotation_final_tr[granulosa_sub$coarse_annotation_final_tr=='Granulosa_sq_transitioning'] = 'Granulosa_sq_transitioning_n'
table(granulosa_sub$coarse_annotation_final_tr,granulosa_sub$coarse_annotation_archr_clean_old)
granulosa_sub$coarse_annotation_final_tr[granulosa_sub$coarse_annotation_archr_clean_old=='Granulosa_sq_transitioning'] = 'Granulosa_sq_transitioning_o'
par(mar=c(0,0,1,20))
plotVisium(granulosa_sub@embeddings$TileMatrix_LSI_harmony_UMAP$df,granulosa_sub$coarse_annotation_final_tr,t='xy',show.cluster.sizes = T)


markerPeaks_tr = getMarkerFeatures(
  ArchRProj = granulosa_sub, 
  useMatrix = "PeakMatrix", 
  groupBy = "coarse_annotation_final_tr",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = 'Granulosa_sq',
  useGroups = c('Granulosa_sq_transitioning_o','Granulosa_sq_transitioning_n'),
  testMethod = "wilcoxon"
)

sgn = assay(markerPeaks,'FDR')   <0.1 & abs(assay(markerPeaks,'Log2FC'))   >0.5
sgn = assay(markerPeaks_tr,'FDR')<0.1 & abs(assay(markerPeaks_tr,'Log2FC'))>0.5
apply(sgn,2,sum)

# saveRDS(markerPeaks,paste0(granulosa_sub@projectMetadata$outputDirectory,'/markerPeaks.rds'))
markerPeaks = readRDS(paste0(granulosa_sub@projectMetadata$outputDirectory,'/markerPeaks.rds'))
peaks_with_genes_ = read.csv('work/pycistopic/results_pycistopic_call_peaks/consensus_peaks_with_genes.csv',row.names = 1)

peaks = rowData(markerPeaks)
mean = assay(markerPeaks,'Mean')
fdr = assay(markerPeaks,'FDR')
l2fc = assay(markerPeaks,'Log2FC')
pv = assay(markerPeaks,'Pval')

rownames(pv) = rownames(mean) = rownames(fdr) = rownames(l2fc) = paste0(peaks$seqnames,':',peaks$start,'-',peaks$end) 
peaks_with_genes = granulosa_sub@peakSet
names(peaks_with_genes) = NULL
peaks_with_genes = as.data.frame(peaks_with_genes)
rownames(peaks_with_genes) = paste0(peaks_with_genes$seqnames,':',peaks_with_genes$start,'-',peaks_with_genes$end)

peaks_with_genes = peaks_with_genes[rownames(fdr),]


colnames(fdr) = paste0('fdr_',colnames(fdr))
colnames(l2fc) = paste0('l2fc_',colnames(l2fc))
colnames(pv) = paste0('pv_',colnames(pv))
dap_sub = cbind(peaks_with_genes,fdr,l2fc,pv)
amh=dap_sub[!is.na(dap_sub$nearestGene) & dap_sub$nearestGene=='AMH',]
amh[order(abs(amh$distToTSS)),-5:-4]
dap_sub['chr19:2248841-2249341',]
# write.csv(dap_sub,paste0(granulosa_sub@projectMetadata$outputDirectory,'/markerPeaks.csv'))




# _check couple of examples ############
# on full dataset to have higher coverage
granulosa = loadArchRProject('work/archr2/04_per_class_clean/Granulosa/')

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

amh_peak = 'chr19:2248841-2249341'
sapply(markers_sgn,function(x)amh_peak %in% rownames(x)[1:5])


for(ct in names(markers_sgn)){
  sgn = markers_sgn[[ct]]
  peaks = c(rownames(sgn[sgn[,paste0('l2fc_',ct)]>0,])[1:5],rownames(sgn[sgn[,paste0('l2fc_',ct)]<0,])[1:5])
  #peak = amh_peak
  pdf(paste0(granulosa_sub@projectMetadata$outputDirectory,'/figures/marker_peaks_examples/',ct,'.pdf'),w=10,h=5)
  for(peak in peaks){
    p <- plotBrowserTrack(
      ArchRProj = granulosa, 
      groupBy = "coarse_annotation_final", 
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

# __plot ATAC near AMH ##############
#peak = amh_peak
pdf(paste0(granulosa_sub@projectMetadata$outputDirectory,'/figures/marker_peaks_examples/AMH.pdf'),w=10,h=5)
peak = amh_peak
p <- plotBrowserTrack(
    ArchRProj = granulosa, 
    groupBy = "coarse_annotation_final", 
    region = toGR(peak,2e4),
    #geneSymbol = c("AMH"),
    useGroups = c('Granulosa_sq','Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml'),
    features =  getMarkers(markerPeaks, cutOff = "abs(Log2FC)>1 & FDR<0.05", returnGR = TRUE)["Granulosa_sq_transitioning"],
    normMethod = 'ReadsInTSS',
    tileSize=50,
  )
  
  grid::grid.newpage()
  grid::grid.draw(p)

dev.off()

# _cluster peaks ########
sgn = apply(fdr<0.05 & abs(l2fc)>2,1,sum)>1

p = rowRanges(mtx)
mtx_ = mtx_[rownames(fdr),]
pmtx_ct = as.matrix(visutils::calcColSums(mtx_,mtx$coarse_annotation_final))
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


# _motif enrichment #####
markerPeaks = readRDS(paste0(granulosa_sub@projectMetadata$outputDirectory,'/markerPeaks.rds'))
# granulosa_sub = addMotifAnnotations(ArchRProj = granulosa_sub, motifSet = "cisbp", annoName = "Motif",force = TRUE)
# granulosa_sub = addMotifAnnotations(ArchRProj = granulosa_sub, motifSet = "JASPAR2020", annoName =  "Motif_JASPAR2020",force = TRUE)
getAvailableMatrices(granulosa_sub)

matchesc = getMatches(ArchRProj = granulosa_sub, name = "Motif")
matchesj = getMatches(ArchRProj = granulosa_sub, name = "Motif_JASPAR2020")
rownames(matchesc) = paste0(seqnames(matchesc),':', start(matchesc),'-', end(matchesc))
matchesc = matchesc[rownames(markerPeaks),]

rownames(matchesj) = paste0(seqnames(matchesj),':', start(matchesj),'-', end(matchesj))
matchesj = matchesj[rownames(markerPeaks),]
# saveRDS(matchesj,paste0(granulosa_sub@projectMetadata$outputDirectory,'/matches_jaspar2020.rds'))
# saveRDS(matchesc,paste0(granulosa_sub@projectMetadata$outputDirectory,'/matches_cisbp.rds'))

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
sgn = assay(markerPeaks,'FDR')[,granulosa_celltypes[-1]] < 0.1 &  assay(markerPeaks,'Log2FC')[,granulosa_celltypes[-1]] <= -0.5

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
cpm =  as.matrix(visutils::calcColSums(cpm,mtx$coarse_annotation_final,mean = TRUE)[,colnames(markerPeaks)])

means = as.matrix(visutils::calcColSums(assay(mtx,'PeakMatrix'),mtx$coarse_annotation_final))[,colnames(markerPeaks)]
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
  pv  = assay(tests[[n]],'Pval')
  or  = assay(tests[[n]],'Enrichment')
  
  colnames(fdr) = paste0('fdr_',colnames(fdr))
  colnames(pv)  = paste0('pv_',colnames(pv))
  colnames(or)  = paste0('oddsRatio_',colnames(or))
  
  write.csv(cbind(pv,or,fdr),paste0(granulosa_sub@projectMetadata$outputDirectory,'/enriched_motifs_',n,'_v2.csv'))
}
saveRDS(tests,paste0(granulosa_sub@projectMetadata$outputDirectory,'/enriched_motifs_v2.rds'))

# compare with previous results
tests_old = readRDS('work/archr/03_granulosa_clean_v2_sub/enriched_motifs_v2.rds')
assayNames(tests_old$cisbp_ind_bkg) = assayNames(tests_old$cisbp_cmn_bkg) = assayNames(tests$cisbp_ind_bkg)

ind     = as.data.frame(myGetMarkers(tests$cisbp_ind_bkg,cutOff = "FDR<0.05 & Log2FC > 0.5",unique = F))
ind_old = as.data.frame(myGetMarkers(tests_old$cisbp_ind_bkg,cutOff = "FDR<0.05 & Log2FC > 0.5",unique = F))

o = names(test_sets)
par(mar=c(20,20,1,1))
imageWithText(compareSets(split(ind$name,ind$comparison),
                          split(ind_old$name,ind_old$comparison),fun = 'o')[o,o])
imageWithText(compareSets(split(ind_old$name,ind_old$comparison),
                          split(ind_old$name,ind_old$comparison),fun = 'o')[o,o])

imageWithText(compareSets(split(ind$name,ind$comparison),
                          split(ind$name,ind$comparison),fun = 'o')[o,o])


cmn = intersect(rownames(tests_old$cisbp_ind_bkg),
                rownames(tests$cisbp_ind_bkg))
length(cmn)
imageWithText(cor(assay(tests_old$cisbp_ind_bkg,'Log2FC')[cmn, ],
    assay(tests$cisbp_ind_bkg,'Log2FC')[cmn, ],method = 'p'))

imageWithText(cor(assay(tests_old$cisbp_ind_bkg,'Pval'),method = 's'))

# comp different bkgs
ind = as.data.frame(myGetMarkers(tests$cisbp_ind_bkg,cutOff = "FDR<0.05 & Log2FC > 0.5",unique = F))
cmn = as.data.frame(myGetMarkers(tests$cisbp_cmn_bkg,cutOff = "FDR<0.05 & Log2FC > 0.5",unique = F))


par(mar=c(20,20,1,1))
imageWithText(compareSets(split(ind$name,ind$comparison),
                          split(cmn$name,cmn$comparison),fun = 'j')[o,o])

compareSets(split(ind$name,ind$comparison),
            split(cmn$name,cmn$comparison),fun = function(x,y)length(intersect(x,y)))[o,o]

compareSets(split(ind$name,ind$comparison),
            split(cmn$name,cmn$comparison),fun = function(x,y)length(union(x,y)))[o,o]

# _cPG islands ############
# downloaded from https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=3213068514_wYJAAFFyL2B1d33A6fmP2Xk8rTqr&db=hg38&c=chr7&g=cpgIslandSuper
cpg = read.csv('data/cpgIslandExt.csv.gz')
peaks = rowData(markerPeaks)
peaks = GRanges(seqnames = peaks$seqnames,IRanges(start=peaks$start,end=peaks$end))
cpg_gr = GRanges(seqnames = cpg$chrom,IRanges(start=cpg$chromStart,end=cpg$chromEnd))
ovlps = findOverlaps(peaks,cpg_gr,type='any',select='all',ignore.strand=T)
peaks_in_cpg = 1:length(peaks) %in% ovlps@from

peaks_in_cpg = matrix(peaks_in_cpg,ncol=1)
rownames(peaks_in_cpg) = rownames(markerPeaks)
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
assay(cpg_tests$ind_bkg,'Pval')
assay(cpg_tests$cmn_bkg,'Pval')
assay(cpg_tests$ind_bkg,'BackgroundProporition')
assay(cpg_tests$cmn_bkg,'Enrichment')


i=rbind(bkg_ind=assay(cpg_tests$ind_bkg,'BackgroundProporition'),
        bkg_cmn=assay(cpg_tests$cmn_bkg,'BackgroundProporition'),
        set=assay(cpg_tests$ind_bkg,'CompareProportion'))



pdf(paste0(granulosa_sub@projectMetadata$outputDirectory,'/figures/CpG_in_marker_peaks.pdf'),w=9,h=6)
par(mfrow=c(1,1),las=2,mar=c(10,4,1,1))
names = c('G_tr_down','G_AMG_early_up','G_AMG_ml_up','G_AMG_down_shared')
b=barplot(i,beside = T,names.arg = names,legend.text = c('individ bkg','common bkg','DARs'),ylim=c(0,0.35),ylab='CpG island proportion')
pv = paste0('pv=',format(assay(cpg_tests$ind_bkg,'Pval'),scientific=T,digits=1))
text(b[1,],i[1,]+0.01,pv,adj=c(0,0),xpd=NA,srt=90)

pv = paste0('pv=',format(assay(cpg_tests$cmn_bkg,'Pval'),scientific=T,digits=1))
text(b[2,],i[2,]+0.01,pv,adj=c(0,0),xpd=NA,srt=90)
dev.off()

barplot(c,beside = T,names.arg = names,legend.text = T,main='Common backgrounds (open in at least one G)',ylim=c(0,0.35))

text(colSums(b)/2,apply(c,2,max)+0.01,pv,adj=c(0.5,0),xpd=NA)


lapply(test_sets,function(x)table(x,peaks_in_cpg[,1]))
lapply(bkg_sets,function(x)table(x,peaks_in_cpg[,1]))


# compare old peaks with new #####
markerPeakso = readRDS('work/archr/03_granulosa_clean_v2_sub/markerPeaks.rds')
peaks = rowData(markerPeaks)
peaks_ = GRanges(seqnames = peaks$seqnames,IRanges(start=peaks$start,end=peaks$end))
dim(peakso)

peakso = rowData(markerPeakso)
peakso_ = GRanges(seqnames = peakso$seqnames,IRanges(start=peakso$start,end=peakso$end))
nrow(peaks)

ovlps = findOverlaps(peaks_,peakso_,type='any',select='all',ignore.strand=T,minoverlap = 251)
ovlps = data.frame(n=rownames(peaks)[ovlps@from],
                   o=rownames(peakso)[ovlps@to],
                   ni=ovlps@from,
                   oi=ovlps@to)
dim(ovlps)
all(peaks$seqnames[ovlps$ni]==peakso$seqnames[ovlps$oi])

ovlps$overlap = pmin(peaks$end[ovlps$ni],peakso$end[ovlps$oi]) - pmax(peaks$start[ovlps$ni],peakso$start[ovlps$oi])
hist(ovlps$overlap,100)
table(table(ovlps$ni))

rownames(ovlps) = ovlps$o


markerPeakso_ = markerPeakso[ovlps$o,]
markerPeaks_ = markerPeaks[ovlps$n,]




t = cor(assay(markerPeaks_,'Pval'),
        assay(markerPeakso_,'Pval'),m='s')

imageWithText(t)


rowData(markerPeaks_)
assayNames(markerPeaks_)

mun = myGetMarkers(markerPeaks_,"FDR<0.1 & Log2FC>0.5",unique = F)
muo = myGetMarkers(markerPeakso_,"FDR<0.1 & Log2FC>0.5",unique = F)
muo$name_n = ovlps[muo$name,'n']

imageWithText(compareSets(split(mun$name,mun$comparison),
            split(muo$name_n,muo$comparison)))


# downsample two flovours of Transitioning ##########################
granulosa = loadArchRProject('work/archr2/04_per_class_clean/Granulosa')
granulosa$coarse_annotation_final_tr = granulosa$coarse_annotation_final
granulosa$coarse_annotation_final_tr[granulosa$coarse_annotation_final_tr=='Granulosa_sq_transitioning'] = 'Granulosa_sq_transitioning_n'

table(granulosa$PeakMatrix_LSI_harmony_cl_2,granulosa$coarse_annotation_archr_clean_old)

granulosa$coarse_annotation_final_tr[granulosa$PeakMatrix_LSI_harmony_cl_2=='C15'] = 'Granulosa_sq_transitioning_o'


plotVisium(granulosa@embeddings$TileMatrix_LSI_harmony_UMAP$df,granulosa$coarse_annotation_final_tr,t='xy',label.clusters = T)

downSample = function(ann,size=min(table(ann)),seed=1234){
  set.seed(seed)
  inx = seq_along(ann)
  inx = unlist(lapply(split(inx,ann),sample,size=size))
  inx
}

cells = granulosa@cellColData
cells = cells[cells$coarse_annotation_final_tr!='Granulosa_AMH_antral_HedgLow',]


inx=downSample(cells$coarse_annotation_final_tr)
table(cells$coarse_annotation_final_tr[inx])

granulosa_sub = subsetArchRProject(
  ArchRProj = granulosa,
  cells = rownames(cells)[inx],
  outputDirectory = "work/archr2/05_per_class_clean_sub/trGranulosa",
  dropCells = TRUE,
  force = TRUE
)

#saveArchRProject(granulosa_sub,outputDirectory = "work/archr2/05_per_class_clean_sub/trGranulosa")
#granulosa_sub = loadArchRProject(path="work/archr2/05_per_class_clean_sub/trGranulosa")
table(granulosa_sub$coarse_annotation_final)
table(granulosa_sub$coarse_annotation_final_tr)

# plot AMH expression ######
amh = gsm['AMH',granulosa_sub$cellNames]
ord = c('Granulosa_sq','Granulosa_sq_transitioning_o','Granulosa_sq_transitioning_n','Granulosa_AMH_early','Granulosa_AMH_ml')

pdf(paste0(granulosa_sub@projectMetadata$outputDirectory,'/figures/AMH_expression.pdf'),w=13,h=7)
par(mfrow=c(2,3),mar=c(0,0,1,10),oma=c(1,1,1,1))
plotVisium(granulosa_sub@embeddings$TileMatrix_LSI_harmony_UMAP$df,granulosa_sub$coarse_annotation_final_tr,t='xy',show.cluster.sizes = T)
plotVisium(granulosa_sub@embeddings$TileMatrix_LSI_harmony_UMAP$df,granulosa_sub$nFrags,t='xy',show.cluster.sizes = T,zfun =log,main='nFrags')
plotVisium(granulosa_sub@embeddings$TileMatrix_LSI_harmony_UMAP$df,amh,t='xy',zfun=log1p,cex=0.6,order.points.by.z = T,legend.args = list(title='AMH cpm'))
par(mar=c(13,4,1,0),las=3)
boxplot(split(log1p(amh), granulosa_sub$coarse_annotation_final_tr)[ord],ylab='log1p AMH cpm')
barplot(tapply(amh, granulosa_sub$coarse_annotation_final_tr,mean)[ord],ylab='mean AMH cpm')
dev.off()


# _check peak matrix ##################
getAvailableMatrices(granulosa_sub)
mtx = getMatrixFromProject(granulosa_sub,'PeakMatrix')
mtx_ = unwrapAssay(mtx)

gene.descr[gene.descr$name=='AMH',]
peaks = as.data.frame(rowRanges(mtx))
f = peaks$seqnames=='chr19' & peaks$start>2248000 & peaks$end<2250000
peaks[f,]
table(mtx_[f,])

# chr19:2248841-2249341 vs
# chr19:2248950-2249450
peak_cnts = readRDS('work/pycistopic/results_pycistopic_call_peaks/peak_chr19_2248950-2249450_cov.rds')
cmn = intersect(colnames(mtx),names(peak_cnts))

# looks ok
table(pycnt=peak_cnts[cmn],archr_new = mtx_[f,cmn])


# _marker peaks #############
markerPeaks = getMarkerFeatures(
  ArchRProj = granulosa_sub, 
  useMatrix = "PeakMatrix", 
  groupBy = "coarse_annotation_final_tr",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = 'Granulosa_sq',
  useGroups = c('Granulosa_sq_transitioning_o','Granulosa_sq_transitioning_n','Granulosa_AMH_early','Granulosa_AMH_ml'),
  testMethod = "wilcoxon"
)


markerPeaks_sq = getMarkerFeatures(
  ArchRProj = granulosa_sub, 
  useMatrix = "PeakMatrix", 
  groupBy = "coarse_annotation_final_tr",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = c('Granulosa_sq_transitioning_o','Granulosa_sq_transitioning_n','Granulosa_AMH_early','Granulosa_AMH_ml'),
  useGroups = 'Granulosa_sq',
  testMethod = "wilcoxon"
)

markerPeaks = cbind(markerPeaks_sq,markerPeaks)
dim(markerPeaks)


up = colSums(assay(markerPeaks,'Log2FC') >  1 & assay(markerPeaks,'FDR') < 0.05)
dw = colSums(assay(markerPeaks,'Log2FC') <  -1 & assay(markerPeaks,'FDR') < 0.05)

up4 = colSums(assay(markerPeaks,'Log2FC') >  4 & assay(markerPeaks,'FDR') < 0.05)
dw4 = colSums(assay(markerPeaks,'Log2FC') <  -4 & assay(markerPeaks,'FDR') < 0.05)


pdf(paste0(granulosa_sub@projectMetadata$outputDirectory,'/figures/marker_peaks_cnt.pdf'),w=8,h=5)
par(mfrow=c(1,2),mar=c(14,4,1,1),las=2)
barplotWithText(log10(up+1),up,main='Peaks up',ylab='log10',col='#EEEEEE',border=NA)
barplotWithText(log10(up4+1),up4,ylab='log10',col='#888888',den=20,add=T,border=NA)
barplotWithText(log10(dw),dw,main='Peaks down',ylab='log10',col='#EEEEEE',border=NA)
barplotWithText(log10(dw4),dw4,ylab='log10',col='#888888',den=20,add=T,border=NA)
legend(grconvertX(0,'nfc','user'),grconvertY(0,'nfc','user'),fill=c('#EEEEEE','#888888'),den=c(-1,20),legend = c('l2FC > 1','l2FC > 4'),xpd=NA,yjust = 0)
dev.off()


p = rowData(markerPeaks)
rownames(markerPeaks) = paste0(p$seqnames,':',p$start,'-',p$end)

# _save markers ###########
# saveRDS(markerPeaks,paste0(granulosa_sub@projectMetadata$outputDirectory,'/markerPeaks.rds'))
markerPeaks = readRDS(paste0(granulosa_sub@projectMetadata$outputDirectory,'/markerPeaks.rds'))
peaks_with_genes_ = read.csv('work/pycistopic/results_pycistopic_call_peaks/consensus_peaks_with_genes.csv',row.names = 1)

peaks = rowData(markerPeaks)
mean = assay(markerPeaks,'Mean')
fdr = assay(markerPeaks,'FDR')
l2fc = assay(markerPeaks,'Log2FC')
pv = assay(markerPeaks,'Pval')

rownames(pv) = rownames(mean) = rownames(fdr) = rownames(l2fc) = paste0(peaks$seqnames,':',peaks$start,'-',peaks$end) 
peaks_with_genes = granulosa_sub@peakSet
names(peaks_with_genes) = NULL
peaks_with_genes = as.data.frame(peaks_with_genes)
rownames(peaks_with_genes) = paste0(peaks_with_genes$seqnames,':',peaks_with_genes$start,'-',peaks_with_genes$end)

peaks_with_genes = peaks_with_genes[rownames(fdr),]


colnames(fdr) = paste0('fdr_',colnames(fdr))
colnames(l2fc) = paste0('l2fc_',colnames(l2fc))
colnames(pv) = paste0('pv_',colnames(pv))
dap_sub = cbind(peaks_with_genes,fdr,l2fc,pv)
amh=dap_sub[!is.na(dap_sub$nearestGene) & dap_sub$nearestGene=='AMH',]
amh[order(abs(amh$distToTSS)),-5:-4]
dap_sub['chr19:2248841-2249341',]
# write.csv(dap_sub,paste0(granulosa_sub@projectMetadata$outputDirectory,'/markerPeaks.csv'))

# _motif enrichment #####
markerPeaks = readRDS(paste0(granulosa_sub@projectMetadata$outputDirectory,'/markerPeaks.rds'))
# granulosa_sub = addMotifAnnotations(ArchRProj = granulosa_sub, motifSet = "cisbp", annoName = "Motif",force = TRUE)
# granulosa_sub = addMotifAnnotations(ArchRProj = granulosa_sub, motifSet = "JASPAR2020", annoName =  "Motif_JASPAR2020",force = TRUE)
getAvailableMatrices(granulosa_sub)

matchesc = getMatches(ArchRProj = granulosa_sub, name = "Motif")
matchesj = getMatches(ArchRProj = granulosa_sub, name = "Motif_JASPAR2020")
rownames(matchesc) = paste0(seqnames(matchesc),':', start(matchesc),'-', end(matchesc))
matchesc = matchesc[rownames(markerPeaks),]

rownames(matchesj) = paste0(seqnames(matchesj),':', start(matchesj),'-', end(matchesj))
matchesj = matchesj[rownames(markerPeaks),]
# saveRDS(matchesj,paste0(granulosa_sub@projectMetadata$outputDirectory,'/matches_jaspar2020.rds'))
# saveRDS(matchesc,paste0(granulosa_sub@projectMetadata$outputDirectory,'/matches_cisbp.rds'))

# _motif enchr tests  #############
# __ define sets ################
motif_enrich_tests = list()
test_sets = bkg_sets = list()


# down in Granulosa_sq_transitioning_o
ct ="Granulosa_sq_transitioning_o"
name='Granulosa_sq_transitioning_o_down'
test_sets[[name]] = assay(markerPeaks,'FDR')[,ct] < 0.1 &  assay(markerPeaks,'Log2FC')[,ct] <= -0.5
bkg_sets[[name]] = assay(markerPeaks,'Mean')[,'Granulosa_sq']>0.01

# down in Granulosa_sq_transitioning_o
ct ="Granulosa_sq_transitioning_n"
name='Granulosa_sq_transitioning_n_down'
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
granulosa_celltypes = colnames(markerPeaks)
ct='Granulosa_AMH_ml'
name = 'Granulosa_AMH_ml_down_shared_vs_only'
sgn = assay(markerPeaks,'FDR')[,granulosa_celltypes[-1]] < 0.1 &  assay(markerPeaks,'Log2FC')[,granulosa_celltypes[-1]] <= -0.5

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
cpm =  as.matrix(visutils::calcColSums(cpm,mtx$coarse_annotation_final_tr,mean = TRUE)[,colnames(markerPeaks)])

means = as.matrix(visutils::calcColSums(assay(mtx,'PeakMatrix'),mtx$coarse_annotation_final_tr))[,colnames(markerPeaks)]
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
  pv  = assay(tests[[n]],'Pval')
  or  = assay(tests[[n]],'Enrichment')
  
  colnames(fdr) = paste0('fdr_',colnames(fdr))
  colnames(pv)  = paste0('pv_',colnames(pv))
  colnames(or)  = paste0('oddsRatio_',colnames(or))
  
  write.csv(cbind(pv,or,fdr),paste0(granulosa_sub@projectMetadata$outputDirectory,'/enriched_motifs_',n,'_v2.csv'))
}
saveRDS(tests,paste0(granulosa_sub@projectMetadata$outputDirectory,'/enriched_motifs_v2.rds'))

# _cPG islands ############
# downloaded from https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=3213068514_wYJAAFFyL2B1d33A6fmP2Xk8rTqr&db=hg38&c=chr7&g=cpgIslandSuper
cpg = read.csv('data/cpgIslandExt.csv.gz')
peaks = rowData(markerPeaks)
peaks = GRanges(seqnames = peaks$seqnames,IRanges(start=peaks$start,end=peaks$end))
cpg_gr = GRanges(seqnames = cpg$chrom,IRanges(start=cpg$chromStart,end=cpg$chromEnd))
ovlps = findOverlaps(peaks,cpg_gr,type='any',select='all',ignore.strand=T)
peaks_in_cpg = 1:length(peaks) %in% ovlps@from

peaks_in_cpg = matrix(peaks_in_cpg,ncol=1)
rownames(peaks_in_cpg) = rownames(markerPeaks)
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
assay(cpg_tests$ind_bkg,'Pval')
assay(cpg_tests$cmn_bkg,'Pval')
assay(cpg_tests$ind_bkg,'BackgroundProporition')
assay(cpg_tests$cmn_bkg,'Enrichment')


i=rbind(bkg_ind=assay(cpg_tests$ind_bkg,'BackgroundProporition'),
        bkg_cmn=assay(cpg_tests$cmn_bkg,'BackgroundProporition'),
        set=assay(cpg_tests$ind_bkg,'CompareProportion'))



pdf(paste0(granulosa_sub@projectMetadata$outputDirectory,'/figures/CpG_in_marker_peaks.pdf'),w=9,h=6)
par(mfrow=c(1,1),las=2,mar=c(10,4,1,1))
names = c('G_tr_down_o','G_tr_down_n','G_AMG_early_up','G_AMG_ml_up','G_AMG_down_shared')
b=barplot(i,beside = T,names.arg = names,legend.text = c('individ bkg','common bkg','DARs'),ylim=c(0,0.35),ylab='CpG island proportion')
pv = paste0('pv=',format(assay(cpg_tests$ind_bkg,'Pval'),scientific=T,digits=1))
text(b[1,],i[1,]+0.01,pv,adj=c(0,0),xpd=NA,srt=90)

pv = paste0('pv=',format(assay(cpg_tests$cmn_bkg,'Pval'),scientific=T,digits=1))
text(b[2,],i[2,]+0.01,pv,adj=c(0,0),xpd=NA,srt=90)
dev.off()

barplot(c,beside = T,names.arg = names,legend.text = T,main='Common backgrounds (open in at least one G)',ylim=c(0,0.35))

text(colSums(b)/2,apply(c,2,max)+0.01,pv,adj=c(0.5,0),xpd=NA)


lapply(test_sets,function(x)table(x,peaks_in_cpg[,1]))
lapply(bkg_sets,function(x)table(x,peaks_in_cpg[,1]))



pdf(paste0(granulosa_sub@projectMetadata$outputDirectory,'/figures/Granulosa_tr_old_vs_new.pdf'),w=13,h=12)
par(mfrow=c(2,3),tcl=-0.2,mgp=c(3.3,0.3,0),mar=c(12,5,1,1),oma=c(0,10,0,0),bty='n',las=3)
granulosa_celltypes = colnames(markerPeaks)
imageWithText(table(granulosa_sub$coarse_annotation_final_tr,granulosa_sub$Sample)[granulosa_celltypes,])
boxplot(split(granulosa_sub$nFrags , granulosa_sub$coarse_annotation_final_tr)[granulosa_celltypes],xlab='',ylab='nFrags',main='nFrags',box=F)
boxplot(split(granulosa_sub$TSSEnrichment , granulosa_sub$coarse_annotation_final_tr)[granulosa_celltypes],xlab='',ylab='TSSEnrichment',main='TSSEnrichment',box=F)
boxplot(split(granulosa_sub$BlacklistRatio , granulosa_sub$coarse_annotation_final_tr)[granulosa_celltypes],xlab='',ylab='BlacklistRatio',main='BlacklistRatio',box=F)
par(mar=c(0,0,1,6))
plotVisium(granulosa_sub@embeddings$TileMatrix_LSI_harmony_UMAP$df,
           granulosa_sub$nFrags,t='xy',zfun =identity,main='nFrags',xlab='',ylab='',label.clusters = granulosa_sub$coarse_annotation_final_tr)
plotVisium(granulosa_sub@embeddings$TileMatrix_LSI_harmony_UMAP$df,
           granulosa_sub$nFrags,t='xy',zfun =log,main='nFrags',xlab='',ylab='',label.clusters = granulosa_sub$coarse_annotation_final_tr)
dev.off()


