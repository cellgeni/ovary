paths = .libPaths()
.libPaths(paths[-grep('/nfs/cellgeni/R',paths)])
library(Seurat)

samples = read.csv('actions/samples.csv')

for(s in samples$sample_id){
  print(s)
  m = Matrix::readMM(paste0('data/',s,'/filtered_peak_bc_matrix/matrix.mtx'))
  b = readLines(paste0('data/',s,'/filtered_peak_bc_matrix/barcodes.tsv'))
  r = data.frame(barcode=b,atac_fragments=colSums(m))
  write.csv(r,paste0('data/',s,'/per_barcode_metrics.csv'),row.names = F,quote = F)
}
