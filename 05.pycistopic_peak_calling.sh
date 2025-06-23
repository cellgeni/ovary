
cd /lustre/scratch127/cellgen/cellgeni/tickets/tic-3942/work/pycistopic
nextflow ../../actions/nf-atac/main.nf \
 --sample_table ../../actions/samples.csv \
 --celltypes ../../work/scanvi_out/combined_gene_matrix_plus_ref_clean_scanvi_raw_obs.csv \
 --output_dir results_pycistopic_call_peaks \
 --fragments_filename fragments.tsv.gz \
 --callPeaks \
 --inferConsensus \
 -resume
