cd work
nextflow ../actions/nf-atac/main.nf \
 -entry snapatac2 \
 --sample_table  ../actions/samples.csv \
 --min_counts 3500 \
 --max_counts 10000000 \
 --output_dir results_snapatac2_combine \
 --fragments_filename fragments.tsv.gz \
 -resume
