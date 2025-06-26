require(seqminer)

getCoverageForSample = function(fragment_file,barcodes2celltype,region,dedupl=TRUE){
  frags = seqminer::tabix.read(fragment_file,tabixRange = paste0(region$chr,':',region$start,'-',region$end))
  frags = as.data.frame(do.call(rbind,strsplit(frags,'\t')))
  colnames(frags) = c('chr_id','start','end','barcode','count')
  if(dedupl)
    frags$count = 1
    
  for(i in c(2,3,5))
    frags[[i]] = as.numeric(frags[[i]]) 
  
  frags = frags[frags$barcode %in% names(barcodes2celltype),]
  frags = split(frags,barcodes2celltype[frags$barcode])
  frags = lapply(frags,countCoverage,region=region)
  return(frags)
}

countCoverage = function(frags,region){
  frags$start = pmax(1,frags$start - region$start + 1)
  frags$end = pmin(region$end-region$start+1,frags$end - region$start + 1)
  
  cov = rep(0,region$end-region$start+1)
  for(i in seq_len(nrow(frags))){
    inx = frags$start[i]:frags$end[i]
    cov[inx] = cov[inx] + frags$count[i]
  }
  
  return(list(x=region$start:region$end,
              chr = region$chr,
              cov=cov))
}


getCoverage = function(fragments_paths,barcodes,region,
                       samples=unique(fragments_paths$sample_id),
                       celltypes=unique(barcodes$celltype),
                       dedupl=TRUE,
                       margin=0){
  region$start = region$start-margin
  region$end = region$end+margin
  fragments_paths = fragments_paths[fragments_paths$sample_id %in% samples,]
  rownames(fragments_paths) = fragments_paths$sample_id
  barcodes = barcodes[barcodes$sample_id %in% fragments_paths$sample_id & barcodes$celltype %in% celltypes,]
  
  samples = unique(barcodes$sample_id)
  celltypes=unique(barcodes$celltype)
  covs = list()
  for(sample in samples){
    barcodes2celltype = barcodes[barcodes$sample_id == sample,]
    barcodes2celltype = setNames(barcodes2celltype$celltype,barcodes2celltype$barcode)
    sample_covs = getCoverageForSample(fragments_paths[sample,'fragment_file'],
                         barcodes2celltype,
                         region,
                         dedupl = dedupl)
    for(celltype in names(sample_covs)){
      if(is.null(covs[[celltype]])){
        covs[[celltype]] = sample_covs[[celltype]]
      }else{
        covs[[celltype]]$cov = covs[[celltype]]$cov + sample_covs[[celltype]]$cov
      }
    }
  }
  return(covs)
}

plotCoverage = function(c,fill='gray',border=NA,ylab='#reads',ylim=range(c$cov,0,1),
                        new=TRUE,xlim=range(c$x),...){
  x = c$x
  y = c$cov
  f = x >= xlim[1] & x <= xlim[2]
  x = x[f]
  y = y[f]
  if(new)
    plot(x,y,t='n',ylab=ylab,ylim=ylim,...)
  polygon(c(x[1],x,x[length(x)]),c(0,y,0),col = fill,border = border)
}

plotCoverages = function(covs,
                         norm_factors=setNames(rep(1,length(covs)),names(covs)),
                         gtf,
                         fill='gray',fill_mark='red',border=NA,ylab='#reads',
                         region2mark=NULL,ylim=NULL,...){
  ylim_ = c(0,1e-10)
  for(n in names(covs)){
    covs[[n]]$cov = covs[[n]]$cov/norm_factors[n]
    ylim_ = range(ylim_,covs[[n]]$cov)
  }
  if(is.null(ylim))
    ylim = ylim_
  
  par(mfrow=c(length(covs)+1,1),mar=c(0,4,1,1),bty='n')
  for(n in names(covs)){
    plotCoverage(covs[[n]],fill=fill,border = border,ylab=ylab,main=n,ylim=ylim,...)
    if(!is.null(region2mark))
      plotCoverage(covs[[n]],fill=fill_mark,border = border,ylab=ylab,main=n,
                   ylim=ylim,xlim=c(region2mark$start,region2mark$end),new=FALSE,...)
  }
  par(mar=c(4,4,1,1))
  tr = gtf[gtf$chr_id==covs[[1]]$chr & gtf$start<covs[[1]]$x[length(covs[[1]]$x)] & gtf$stop>covs[[1]]$x[1] , ]
  plotTranscripts(tr,xlim=range(covs[[1]]$x))
  
}

pasreCoors = function(coors){
  coors = strsplit(coors,'[:-]')
  coors = as.data.frame(do.call(rbind,coors))
  colnames(coors) = c('chr','start','end')
  coors$start = as.numeric(coors$start)
  coors$end = as.numeric(coors$end)
  coors
}

plotTranscripts = function (a, ylim = c(0, length(unique(a$transcript_id))), 
                            xlim = c(ifelse(a$strand[1] == "+", min(a$start), max(a$stop)), 
                                     ifelse(a$strand[1] == "+", max(a$stop), min(a$start))), 
                            xlab = a$chr_id[1], new = TRUE, 
                            yspace = 0.8, exon.col = "black", cds.col = "black", text.cex = 1,...) {
  if (!is.na(exon.col)) 
    a$exon.col = exon.col
  if (!is.na(cds.col)) 
    a$cds.col = cds.col
  transc = split(a, a$transcript_id)
  transc = transc[order(sapply(transc, function(x) {
    max(x$stop) - min(x$start)
  }))]
  if (new) 
    plot(1, t = "n", xlim = xlim, ylim = ylim, yaxt = "n", 
         ylab = "", xlab = xlab, ...)
  ystep = (ylim[2] - ylim[1])/length(transc)
  for (i in 1:length(transc)) {
    y = ylim[1] + ystep * i - ystep/2
    t = transc[[i]]
    arrowSegment(min(t$start),y,max(t$stop),y,width=ystep/8,angle = ifelse(t$strand[1]=='+',35,ifelse(t$strand[1]=='-',145,90)))
    f = t$feature == "exon"
    if (sum(f) > 0) 
      rect(t$start[f], y - ystep/2 * yspace, t$stop[f], 
           y + ystep/2 * yspace, col = "white", border = t$exon.col[f])
    f = t$feature == "CDS"
    if (sum(f) > 0) 
      rect(t$start[f], y - ystep/2 * yspace, t$stop[f], 
           y + ystep/2 * yspace, col = t$cds.col[f], border = t$cds.col[f])
  }
  text(par("usr")[1], seq(ylim[1] + ystep/2, by = ystep, length.out = length(transc)), 
       sapply(transc, function(x) x$gene_name[1]), adj = c(1, 
                                                                 0.5), xpd = T, cex = text.cex)
}

gconvertL = function(x,what,...){
  r = do.call(paste0('grconvert',toupper(what)),list(c(0,x),...))
  r[2]-r[1]
}

arrowSegment = function(x0,y0,x1,y1,den=10,angle=45,width=0.3,...){
  segments(x0,y0,x1,y1)#,...)
  
  tick_space = abs(gconvertL(1/den,what='x',from='inch',to='user'))
  tick_pos = seq(x0,x1,by=tick_space)
  tick_shift = gconvertL(gconvertL(width,'y',from='user',to='inch'),'x','inch','user')/tan(angle/180*pi)
  segments(tick_pos,y0,
           x1=tick_pos-tick_shift,
           y0+width)
  
  segments(tick_pos,y0,
           x1=tick_pos-tick_shift,
           y0-width)
}
#select longest protein coding transcript per gene
# gtf = readRDS('/nfs/cellgeni/pasham/code/nf-scsajr/ref/human_2020A_chr/gtf.rds')
# gtf = gtf[gtf$transcript_type=='protein_coding',]
# gtf = split(gtf,gtf$gene_id)
# gtf = lapply(gtf,function(g){
#   t = g[g$feature=='transcript',]
#   tid_longest = t$transcript_id[order(t$stop-t$start,decreasing = T)[1]]
#   rbind(g[g$feature=='gene',],
#         g[g$transcript_id==tid_longest,])
# })
# gtf = do.call(rbind,gtf)
# dim(gtf)
# saveRDS(gtf,'gtf_longest_tr.rds')

# or keep all exons
# gtf = readRDS('/nfs/cellgeni/pasham/code/nf-scsajr/ref/human_2020A_chr/gtf.rds')
# gtf = gtf[gtf$feature %in% c('exon','CDS'),c('chr_id','feature','start','stop','strand','gene_id','gene_name')]
# dim(gtf)
# gtf = unique(gtf)
# gtf$transcript_id = gtf$gene_id
# 
# saveRDS(gtf,'gtf_all_exons_on_one_transc.rds')
# 
# a=gtf[gtf$chr_id=='chr1' & gtf$start>=1e5 & gtf$stop<5e5,]
# plotTranscripts(a)
