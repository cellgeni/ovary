nextflow ~/nfs/projects/2504.ovary_atac/src/nf-atac/main.nf \
 -entry snapatac2 \
 --sample_table /lustre/scratch127/cellgen/cellgeni/pasham/data/2504.ovary_atac/atac/raw/samples.csv \
 --min_counts 3500 \
 --max_counts 10000000 \
 --celltypes /lustre/scratch127/cellgen/cellgeni/pasham/data/2504.ovary_atac/atac/scanvi_out/celltypes_01.csv \
 --output_dir /lustre/scratch127/cellgen/cellgeni/pasham/data/2504.ovary_atac/atac/results2 \
 --remove_doublets True \
 -resume
