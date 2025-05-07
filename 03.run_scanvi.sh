#! /bin/bash -e
#BSUB -G cellgeni
#BSUB -J scanvi
#BSUB -o logs/%J.scanvi.log
#BSUB -e logs/%J.scanvi.err
#BSUB -n 1
#BSUB -M80000
#BSUB -R "span[hosts=1] select[mem>80000] rusage[mem=80000]"

# for gpu-normal
#BSUB -q gpu-normal
#BSUB -gpu "mode=shared:j_exclusive=yes:gmem=32000:num=1"

WDIR=`pwd -P`
module load cellgen/singularity

IMAGE=/nfs/cellgeni/singularity/images/scvi-1.1.2.sif # _metrics.sif # with scvi v1.3.0, but it has no scikit-misc

  
singularity exec --nv --bind /lustre,/nfs $IMAGE /bin/bash -c "nvidia-smi;cd ${WDIR}; \
 ./src/ovary/bin/run_scanvi.py \
  --h5ad_path data.lustre/atac/scanvi_out/combined_gene_matrix_plus_ref.h5ad \
  --h5ad_out data.lustre/atac/scanvi_out/combined_gene_matrix_plus_ref_scanvi_ds-dd_e1000.h5ad \
  --batch_key dataset \
  --categorical_covariate_keys dataset_donor \
  --celltype_key coarse_annotation \
  --n_top_genes 5000 \
  --max_epochs 1000 \
  --lr 0.0005"