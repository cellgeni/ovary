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
        description="Builds mallet model"
    )
    parser.add_argument(
        "--cistopic_pkl_in", 
        type=str,  
        help="Path cistopic obj",
    )
    
    parser.add_argument(
        "--pkl_out", 
        type=str,  
        help="file to save model",
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
    args = parser.parse_args()
    print(args)

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

    pickle.dump(
        models,
        open(args.pkl_out, "wb")
        )

if __name__ == "__main__":
    main()
