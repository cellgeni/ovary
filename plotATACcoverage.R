.libPaths('/software/cellgen/cellgeni/R_ATAC/') # dependencies are installed here
source('/nfs/cellgeni/tickets/tic-3942/actions/ovary/bin/plotATACCoverage.R') 

# define paths to inputs data
path = '/lustre/scratch127/cellgen/cellgeni/tickets/tic-3942/'
# sample_id to fragment file link
fragments_paths = read.csv(paste0(path,'actions/samples.csv'))
fragments_paths$fragment_file = paste0(fragments_paths$filedir,'/fragments.tsv.gz')

# celltype annotation
barcodes = read.csv(paste0(path,'work/scanvi_out/combined_gene_matrix_plus_ref_clean_scanvi_raw_obs.csv'))
# gene coordinates (longest protein coding isoform per gene based on 2020A 10x reference)
# gtf=readRDS(paste0(path,'gtf_longest_tr.rds')) # to plot longest protein coding isoform
gtf=readRDS(paste0(path,'gtf_all_exons_on_one_transc.rds')) # to plot all exons of the gene on one line
# cell*peak counts. loaded only for normalization normalization
atac = schard::h5ad2list(paste0(path,'work/results_snapatac2_call_peaks/subset_adatas/peak_mat_granulosa.h5ad')) # granulosa only (faster to load than whole object)
# atac = schard::h5ad2list(paste0(path,'work/results_snapatac2_call_peaks/subset_adatas/peak_mat.h5ad')) # all celltypes
celltype_totals = sapply(split(atac$obs$total_counts,atac$obs$celltype),sum)


# peak to plot
peak = pasreCoors('chr19:2240846-2241863')
# celltypes to plots
cts = c('Granulosa_sq','Granulosa_sq_atr','Granulosa_sq_transitioning','Granulosa_AMH_early','Granulosa_AMH_ml')
# margin width to include into plot
margin = 5e4
c1=getCoverage(fragments_paths,barcodes,peak,dedupl=FALSE,margin=margin)

plotCoverages(c1[cts],celltype_totals[cts]/1e6,ylab='CPM',xlim=c(2240000,2250000),
              region2mark = peak,
              gtf=gtf,
              xaxt='n',xlab='',ylim=NULL)
