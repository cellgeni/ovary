#!/usr/bin/env python3

# singularity run -B /nfs,/lustre /nfs/cellgeni/singularity/images/scenicplus-fa55dae.sif python
import pycisTopic
import dill
import scanpy as sc
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
import numpy as np

from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)

from pycisTopic.gene_annotation import (
    get_chrom_sizes_and_alias_mapping_from_ucsc
)

import pandas as pd
import tempfile
import os
import pyranges as pr
from pycisTopic.gene_activity import get_gene_activity
from pycisTopic.lda_models import run_cgs_models_mallet
import argparse
import pickle

def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Generates gene score matrix using pycistopic"
    )
    parser.add_argument(
        "--cistopic_pkl_in", 
        type=str,  
        help="Path cistopic obj",
    )
    
    parser.add_argument(
        "--cistopic_pkl_out", 
        type=str,  
        help="Output pkl",
    )
    parser.add_argument(
        '--n_topics', 
        nargs='+', 
        type=int, 
        help='List of number of topics')
    
    parser.add_argument(
        "--ncpu", 
        type=int,  
        help="ncpu",
    )
    
    parser.add_argument(
        "--mem", 
        type=int,  
        help="RAM in G",
    )
    return parser


def main():
    # parse script arguments
    parser = init_parser()
    #args = parser.parse_args()
    args = parser.parse_args(["--cistopic_pkl_in","/lustre/scratch127/cellgen/cellgeni/tickets/tic-3942/work/pycistopic/results_pycistopic_combine/combined_cistopic_object.pkl",
                              "--cistopic_pkl_out","/lustre/scratch127/cellgen/cellgeni/tickets/tic-3942/work/pycistopic/results_pycistopic_combine/combined_cistopic_object_mod.pkl",
                              "--n_topics",'40',
                              '--ncpu','16',
                              '--mem','200'])
    
    with open(args.cistopic_pkl_in, "rb") as f:
        cistopic = dill.load(f)
    
    os.environ['MALLET_MEMORY'] = str(args.mem)+'G'

    # Configure path Mallet
    mallet_path="/opt/Mallet/bin/mallet"
    tmp_path = tempfile.mkdtemp()#'/tmp'
    # Run models
    models=run_cgs_models_mallet(
        cistopic,
        n_topics=args.n_topics,
        n_cpu=args.ncpu,
        n_iter=500,
        random_state=555,
        alpha=50,
        alpha_by_topic=True,
        eta=0.1,
        eta_by_topic=False,
        tmp_path=tmp_path,
        save_path=tmp_path,
        mallet_path=mallet_path,
    )
    cistopic.add_LDA_model(models[0])

    pickle.dump(
        cistopic,
        open(args.cistopic_pkl_out, "wb")
    )

    imputed_acc_obj = impute_accessibility(
        cistopic,
        selected_cells=None,
        selected_regions=None,
        scale_factor=10**6
    )

    pickle.dump(
        imputed_acc_obj,
        open(args.cistopic_pkl_out+"imputed_acc_obj.pkl", "wb")
    )

    #with open(args.cistopic_pkl_out+"imputed_acc_obj.pkl", "rb") as f:
    #    imputed_acc_obj = dill.load(f)
    
    #del cistopic

    chromsizes = get_chrom_sizes_and_alias_mapping_from_ucsc(
       ucsc_assembly="hg38",
       chrom_sizes_and_alias_tsv_filename="hg38.chrom_sizes_and_alias.tsv",
    )

    # it is essential to reload it, as get_chrom_sizes_and_alias_mapping_from_ucsc return polaris that is not compatable with get_gene_activity
    chromsizes = pd.read_table("hg38.chrom_sizes_and_alias.tsv")
    chromsizes.rename({"# ucsc": "Chromosome", "length": "End"}, axis = 1, inplace = True)
    chromsizes["Start"] = 0
    chromsizes = pr.PyRanges(chromsizes[["Chromosome", "Start", "End"]])

    pr_annotation = pd.read_table(
        "actions/nf-atac/reference/hg38_pycistopic_tss.bed"
    ).rename({"Name": "Gene", "# Chromosome": "Chromosome"}, axis = 1)
    pr_annotation["Transcription_Start_Site"] = pr_annotation["Start"]
    pr_annotation = pr.PyRanges(pr_annotation)

    # at least 350G
    gene_act, weigths = get_gene_activity(
        imputed_acc_obj,
        pr_annotation,
        chromsizes,
        use_gene_boundaries=True, # Whether to use the whole search space or stop when encountering another gene
        upstream=[1000, 100000], # Search space upstream. The minimum means that even if there is a gene right next to it
                             # these bp will be taken (1kbp here)
        downstream=[1000,100000], # Search space downstream
        distance_weight=True, # Whether to add a distance weight (an exponential function, the weight will decrease with distance)
        decay_rate=1, # Exponent for the distance exponential funciton (the higher the faster will be the decrease)
        extend_gene_body_upstream=10000, # Number of bp upstream immune to the distance weight (their value will be maximum for
                          #this weight)
        extend_gene_body_downstream=500, # Number of bp downstream immune to the distance weight
        gene_size_weight=False, # Whether to add a weights based on the length of the gene
        gene_size_scale_factor='median', # Dividend to calculate the gene size weigth. Default is the median value of all genes
                          #in the genome
        remove_promoters=False, # Whether to remove promoters when computing gene activity scores
        average_scores=True, # Whether to divide by the total number of region assigned to a gene when calculating the gene
                          #activity score
        scale_factor=1, # Value to multiply for the final gene activity matrix
        extend_tss=[10,10], # Space to consider a promoter
        gini_weight = True, # Whether to add a gini index weigth. The more unique the region is, the higher this weight will be
        return_weights= True, # Whether to return the final weights
        project='Gene_activity') # Project name for the gene activity object
    
    
    gene_act_adata = sc.AnnData(gene_act.mtx.T, obs = cistopic.cell_data.copy(), var = pd.DataFrame(index = gene_act.feature_names))

    
    gene_act_adata.obs_names = gene_act_adata.obs.sample_id + ":"+gene_act_adata.obs.barcode
    gene_act_adata.obs.cisTopic_nr_frag     = gene_act_adata.obs.cisTopic_nr_frag.astype(int)
    gene_act_adata.obs.cisTopic_log_nr_frag = gene_act_adata.obs.cisTopic_log_nr_frag.astype(float)
    gene_act_adata.obs.cisTopic_nr_acc      = gene_act_adata.obs.cisTopic_nr_acc.astype(int)
    gene_act_adata.obs.cisTopic_log_nr_acc  = gene_act_adata.obs.cisTopic_log_nr_acc.astype(float)

    gene_act_adata.write_h5ad(args.cistopic_pkl_out.replace('.pkl','')+"_gene_scores.h5ad")

    run_umap(
        cistopic,
        target  = 'cell', scale=True)
    
    find_clusters(
       cistopic,
        target  = 'cell',
        k = 10,
        res = [0.6, 1.2, 3],
        prefix = 'pycisTopic_',
        scale = True,
        split_pattern = '-'
    )

    pickle.dump(
        cistopic,
        open(args.cistopic_pkl_out, "wb")
    )


        
if __name__ == "__main__":
    main()


