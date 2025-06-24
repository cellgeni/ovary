nextflow ../actions/nf-atac/main.nf \
 -entry snapatac2 \
 --sample_table  ../actions/samples.csv \
 --min_counts 3500 \
 --max_counts 10000000 \
 --celltypes ../work/scanvi_out/combined_gene_matrix_plus_ref_clean_scanvi_raw_obs.csv \
 --output_dir results_snapatac2_call_peaks_replicates \
 --fragments_filename fragments.tsv.gz \
 --remove_doublets True \
 -resume
