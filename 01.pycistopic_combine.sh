# generate dummy celltype annotation (just one celltype)
# groups all cells together as pysictopyc need peaks (and then topics and imputed accesabilities) to calcualte gene scores
echo 'sample_id,barcode,celltype' > work/atac_celltypes_dummy.csv
for s in `cut -f 1 -d ' ' actions/samples_irods.txt`
do
 for b in `cut -f 1 -d ',' data/$s/per_barcode_metrics.csv | grep -v 'barcode'`
 do
	 echo "${s},${b},all"
 done
done >> work/atac_celltypes_dummy.csv

# prepare barcode list 
singularity run --bind /nfs,/lustre /nfs/cellgeni/singularity/images/snareseq_star2.7.10b_bwamem2v2.2.1_archr1.0.2_seurat4.2.1_visutils.sif "Rscript actions/00.pycistopic_prepare_barcode_metrics.R"

cd work/pycistopic
# just generate cistopic object, all cells as one celltype 
nextflow ../../actions/nf-atac/main.nf \
 --sample_table ../../actions/samples.csv \
 --celltypes  ../../work/atac_celltypes_dummy.csv \
 --output_dir results_pycistopic_combine \
 --fragments_filename fragments.tsv.gz \
 --callPeaks \
 --inferConsensus \
 -resume
