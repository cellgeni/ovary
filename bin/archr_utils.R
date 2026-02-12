plotExpAsDots = function(m,f,gids,scale_exp_mar=1,max.cex=1,...){
  m = m[gids,]
  e = as.matrix(visutils::calcColSums(m,f,mean=TRUE))
  n = as.matrix(visutils::calcColSums(m>0,f,mean=TRUE))
  if(!is.na(scale_exp_mar)){
    e = sweep(e,scale_exp_mar,apply(e,scale_exp_mar,mean),'-')
    e = sweep(e,scale_exp_mar,apply(e,scale_exp_mar,sd),'/')
  }
  visutils::dotPlot(n,e,max.cex = max.cex,legend.cex.at = 0:4/4,legend.col.at = NULL,
                    legend.cex.title = 'fraction cells',legend.col.title = 'mean cpm',
                    colfun = function(x) num2col(x, c("yellow", "violet", "black")),...)
}

# modified version of ArchR::getMarkers that returns all assays
myGetMarkers = function (seMarker = NULL, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", 
          n = NULL, returnGR = FALSE,simplify=TRUE,unique=TRUE) {
  assayNames <- names(SummarizedExperiment::assays(seMarker))
  for (an in assayNames) {
    eval(parse(text = paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", 
                             an, "']]")))
  }
  passMat <- eval(parse(text = cutOff))
  for (an in assayNames) {
    eval(parse(text = paste0("rm(", an, ")")))
  }
  if (returnGR) {
    if (S4Vectors::metadata(seMarker)$Params$useMatrix != 
        "PeakMatrix") {
      stop("Only markers can be returned as GRanges when PeakMatrix!")
    }
    rr <- GRanges(rowData(seMarker)$seqnames, IRanges(rowData(seMarker)$start, 
                                                      rowData(seMarker)$end))
    grL <- lapply(seq_len(ncol(passMat)), function(x) {
      idx <- which(passMat[, x])
      rrx <- rr[idx]
      if(nrow(rrx)==0)
        return(NULL)
      rrx[,'comparison'] = rep(colnames(passMat)[x],nrow(rrx))
      for (an in assayNames) {
        rrx[,an] <- SummarizedExperiment::assays(seMarker[idx,])[[an]][, x]
      }
      rrx <- rrx[order(rrx$Pval), , drop = FALSE]
      if (!is.null(n)) {
        if (n < nrow(rrx)) {
          rrx <- rrx[seq_len(n), , drop = FALSE]
        }
      }
      rrx
    }) %>% SimpleList
    names(grL) <- colnames(seMarker)
    grL <- grL[gtools::mixedsort(names(grL))]
    return(grL)
  }
  else {
    markerList <- lapply(seq_len(ncol(passMat)), function(x) {
      idx <- which(passMat[, x])
      rrx <- SummarizedExperiment::rowData(seMarker[idx,])
      rrx[,'comparison'] = rep(colnames(passMat)[x],nrow(rrx))
      rrx[,'name'] = rownames(rrx)
      for (an in assayNames) {
        rrx[,an] <- SummarizedExperiment::assays(seMarker[idx,])[[an]][, x]
      }
      rrx <- rrx[order(rrx$Pval), , drop = FALSE]
      rrx
    })
    
    markerList = do.call(rbind,markerList)
    if(unique){
      markerList = do.call(rbind,lapply(split(markerList,markerList$name),function(x){x[order(x$Log2FC,decreasing = T)[1],]}))
      rownames(markerList) = markerList$name
    }
    markerList = split(markerList,markerList$comparison)
    markerList = lapply(markerList,function(x)x[order(x$Pval),])
    if(!is.null(n)){
      markerList = lapply(markerList,function(x)x[seq_len(min(n,nrow(x))),])
    }
    if(simplify)
      markerList = do.call(rbind,markerList)

    return(markerList)
  }
}

# modified version of peakAnnoEnrichment that can work with custom bg
myPeakAnnoEnrichment = function(p2m,set,bkg,name='set'){
  if(is.null(rownames(p2m)))
    rownames(p2m) = 1:nrow(p2m)
  set = rownames(p2m[set,,drop=FALSE])
  bkg = union(set,rownames(p2m[bkg,,drop=FALSE]))
  if(is(p2m,"SummarizedExperiment"))
    p2m = assay(p2m,'matches')
  p2m = p2m[bkg,,drop=FALSE]
  
  set = bkg %in% set
  pv = l2or = BackgroundFrequency = nBackground = nCompare = CompareFrequency = matrix(0,ncol=1,nrow=ncol(p2m),dimnames = list(colnames(p2m),name))
  nBackground[,] = nrow(p2m)
  nCompare[,] = sum(set)
  BackgroundFrequency[,1] = colSums(p2m)
  CompareFrequency[,1] = colSums(p2m[set,,drop=FALSE])
  for(i in 1:ncol(p2m)){
    pv[i,1] = phyper(CompareFrequency[i,1] - 1, # Number of Successes the -1 is due to cdf integration
                     BackgroundFrequency[i,1], # Number of all successes in background
                     nBackground[i,1] - BackgroundFrequency[i,1], # Number of non successes in background
                     nCompare[i,1], # Number that were drawn
                     lower.tail = FALSE, log.p = FALSE)
  }
  BackgroundProporition = BackgroundFrequency/nBackground
  CompareProportion = CompareFrequency/nCompare
  Enrichment = CompareProportion/BackgroundProporition
  FDR = pv
  FDR[] = apply(pv,2,p.adjust,m='BH')
  res = SummarizedExperiment(
    assays = list(Pval = pv,
                  FDR = FDR,
                  mlog10Padj = -log10(FDR),
                  mlog10p = -log10(pv),
                  Enrichment = Enrichment,
                  Log2FC = log2(Enrichment),
                  BackgroundProporition = BackgroundProporition,
                  nBackground = nBackground,
                  BackgroundFrequency = BackgroundFrequency,
                  CompareProportion = CompareProportion,
                  nCompare = nCompare,
                  CompareFrequency = CompareFrequency
    )
  )
}


myCreateGeneAnnotation = function(genome,gtf_file){
  addArchRGenome(genome)
  gtf = rtracklayer::import(gtf_file)
  
  gtf = gtf[gtf@seqnames %in% levels(getChromSizes()@seqnames),]
  ann = createGeneAnnotation(genome = genome,
                             genes = gtf[gtf$type=='gene',],
                             exons = gtf[gtf$type=='exon',])
  ann$genes$symbol = ann$genes$gene_name
  ann
}


run_rdim_cl = function(p,batch,mtx_names,cl_ress=c(1,2,4,6,10)){
  
  for(m in mtx_names){
    p = addIterativeLSI(
      ArchRProj = p,
      useMatrix = m, 
      name = paste0(m,"_LSI"), 
      iterations = 4, 
      clusterParams = list(resolution = c(2), sampleCells = 10000, maxClusters = 6, n.start = 10), 
      varFeatures = 25000, 
      dimsToUse = 1:30, 
      LSIMethod = 2,
      force = TRUE
    )
    
    p = addHarmony(
      ArchRProj = p,
      reducedDims = paste0(m,"_LSI"),
      name = paste0(m,"_LSI_harmony"),
      groupBy = batch,
      force = TRUE
    )
  }
  
  
  for(dim in names(p@reducedDims)){
    print(dim)
    for(r in cl_ress){
      print(r)
      p <- addClusters(
        input = p,
        reducedDims = dim,
        method = "Seurat",
        name = paste0(dim,'_cl_',r),
        resolution = r, 
        #maxClusters = 25, 
        knnAssign = 20, 
        force = TRUE
      )
    }
    
    p <- addUMAP(
      ArchRProj = p,
      reducedDims = dim,
      name = paste0(dim,'_UMAP'),
      nNeighbors = 30,
      minDist = 0.5,
      metric = "cosine",
      force = TRUE
    )
  }
  p
}

unwrapAssay = function(e,a=SummarizedExperiment::assayNames(e)[1]){
  r = SummarizedExperiment::assay(e,a)
  if('name' %in% colnames(rowData(e))){
    rownames(r) = rowData(e)[,'name']
  }else if('RangedSummarizedExperiment' %in% class(e)){
    rr = as.data.frame(rowRanges(e))
    rownames(r) = paste0(rr$seqnames,':',rr$start,'-',rr$end)
  }
  r
}


compareMtxs = function(a,b){
  rcmn = intersect(rownames(a),rownames(b))
  ccmn = intersect(colnames(a),colnames(b))
  # rownames
  if(nrow(a) == nrow(b) && nrow(a) == length(rcmn))
    print(paste0("rows (a,b,a&b): SAME"))
  else  
    print(paste0("rows (a,b,a&b):", nrow(a),',',nrow(b),',',length(rcmn)))
  # colnames
  if(ncol(a) == ncol(b) && ncol(a) == length(ccmn))
    print(paste0("cols (a,b,a&b): SAME"))
  else  
    print(paste0("cols (a,b,a&b):", ncol(a),',',ncol(b),',',length(ccmn)))
  
  if(length(rcmn)==0 | length(ccmn) ==0){
    print("no common obs")
    return(0)
  }
    
  a = a[rcmn,ccmn]
  b = b[rcmn,ccmn]
  # rowsums
  print(paste0("row sums are same: ",all(rowSums(a)==rowSums(b))))
  # colsums
  print(paste0("col sums are same: ",all(colSums(a)==colSums(b))))
  # identical
  print(paste0("identical:", all(a@i==b@i) & all(a@p==b@p) & all(a@x==b@x)))
}

getMajorCelltype = function(cl,celtype){
  t = table(celtype,cl)
  cl2class = setNames(rownames(t)[apply(t,2,which.max)],colnames(t))
  cl2class[cl]
}


