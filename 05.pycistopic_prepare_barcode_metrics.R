library(Seurat)
samples = read.csv('data.lustre/atac/raw/samples.csv')
for(s in samples$sample_id){
  print(s)
  m = Matrix::readMM(paste0('data.lustre/atac/raw/',s,'/filtered_peak_bc_matrix/matrix.mtx'))
  b = readLines(paste0('data.lustre/atac/raw/',s,'/filtered_peak_bc_matrix/barcodes.tsv'))
  #r = read.csv(paste0('data.lustre/atac/raw/',s,'/per_barcode_metrics.csv'))
  #colnames(r)[2] = 'atac_fragments'
  r = data.frame(barcode=b,atac_fragments=colSums(m))
  write.csv(r,paste0('data.lustre/atac/raw/',s,'/per_barcode_metrics.csv'),row.names = F,quote = F)
}
