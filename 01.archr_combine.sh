# cd work/archr
nextflow ../../actions/nf-atac/main.nf \
 -entry archr \
 --sample_table  ../../actions/samples.csv \
 --min_counts 3500 \
 --output_dir 01_arrows \
 --fragments_filename fragments.tsv.gz \
 --genome hg38 \
 -resume
