
#nextflow run cellgeni/nf-atac -r 25-056 \

nextflow ~/nfs/projects/2504.ovary_atac/src/nf-atac/main.nf \
 --sample_table /lustre/scratch127/cellgen/cellgeni/pasham/data/2504.ovary_atac/atac/raw/samples_pycistopic.csv \
 --celltypes /lustre/scratch127/cellgen/cellgeni/pasham/data/2504.ovary_atac/atac/scanvi_out/celltypes_noslash_01.csv \
 --output_dir pycistopic \
 --fragments_filename fragments.tsv.gz \
 --callPeaks \
 --inferConsensus \
 -resume
