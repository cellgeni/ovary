#!/usr/bin/env python3

import warnings
warnings.filterwarnings("ignore")

import anndata as ad
import pandas as pd
import scanpy as sc
import scvi
import numpy as np
import argparse
import torch
import os

scvi.settings.seed = 0




def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description=""
    )
    parser.add_argument(
        "--h5ad_path", 
        type=str, 
        help="Path to object",
    )
    
    parser.add_argument(
        "--h5ad_out", 
        type=str, 
        help="Path to write output object",
    )
    
    
    parser.add_argument(
        "--batch_key",
        type=str,
        help="Name of batch column in adata.obs"
    )
    
    parser.add_argument(
        "--categorical_covariate_keys",
        type=str,
        help="Name of batch column in adata.obs",
        action='append'
    )
    
    parser.add_argument(
        "--celltype_key",
        type=str,
        help="Name of celltype column in adata.obs"
    )
    
    parser.add_argument(
        "--n_top_genes",
        type=int,
        help="Number of variable genes to use"
    )
    
    parser.add_argument(
        "--max_epochs",
        type=int,
        help="Max epochs",
        default = 1000
    )

    parser.add_argument(
        "--lr",
        type=float,
        help="learning rate",
        default = 0.001
    )
    
    parser.add_argument(
        "--hgv_flavor",
        type=str,
        help="flavor for highly_variable_genes. Data will be log1p-ed for HVG (but not for integration) if not seurat_v3/seurat_v3_paper.",
        default = 'seurat_v3'
    )
    
    return parser



def main():
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()
    
    data = sc.read_h5ad(args.h5ad_path)
    
    sc.pp.filter_genes(data, min_cells=5)
    
    data.layers['original'] = data.X.copy()
    
    # log data for HGV in case flavor is not seurat_v3
    if args.hgv_flavor not in {'seurat_v3','seurat_v3_paper'}:
        sc.pp.log1p(data)
    
    sc.pp.highly_variable_genes(
      data,
      n_top_genes = args.n_top_genes,
      flavor=args.hgv_flavor,
      batch_key=args.batch_key,
      subset=True
    )
    
    # replace original expression back
    data.X = data.layers['original']
    del data.layers['original']
    
    pred_key = args.celltype_key + "_scanvi"
    
    data.obs[pred_key] = 'Unknown'
    ref_idx = ~data.obs[args.celltype_key].isna()
    data.obs[pred_key][ref_idx] = data.obs[args.celltype_key][ref_idx]
    data.obs[pred_key] = data.obs[pred_key].astype('category')
    
    scvi.model.SCVI.setup_anndata(data,
      batch_key=args.batch_key,
      categorical_covariate_keys=args.categorical_covariate_keys)
    
    vae = scvi.model.SCVI(
      data,
      n_layers=2,
      n_latent=30,
      gene_likelihood="nb",
      dispersion="gene-batch"
    )

    if os.path.isdir(args.h5ad_out+"_scvi_mod"):
        vae = scvi.model.SCVI.load(args.h5ad_out+"_scvi_mod",data)
    else:
        vae.train(max_epochs=args.max_epochs, early_stopping=True,plan_kwargs={"lr": args.lr})    
        vae.save(dir_path=args.h5ad_out+"_scvi_mod")
    
    # plot scvi train
    train_elbo = vae.history["elbo_train"][1:]
    test_elbo = vae.history["elbo_validation"]
    ax = train_elbo.plot()
    test_elbo.plot(ax=ax)
    ax.figure.savefig(args.h5ad_out+"_scvi_mod/training.pdf")
    
    data.obsm["X_scVI"] = vae.get_latent_representation(data)
    sc.pp.neighbors(data, use_rep="X_scVI")
    sc.tl.umap(data)
    data.obsm['X_scVI_umap'] = data.obsm['X_umap']

    # scANVI
    lvae = scvi.model.SCANVI.from_scvi_model(
      vae,
      adata=data,
      labels_key=pred_key,
      unlabeled_category="Unknown",
      #var_activation=torch.nn.functional.softplus,
    )
    lvae.train(max_epochs=args.max_epochs, n_samples_per_label=100,early_stopping=True,plan_kwargs={"lr": args.lr})
    lvae.save(dir_path=args.h5ad_out+"_scanvi_mod")
    
    train_elbo = lvae.history["elbo_train"][1:]
    test_elbo = lvae.history["elbo_validation"]
    ax = train_elbo.plot()
    test_elbo.plot(ax=ax)
    ax.figure.savefig(args.h5ad_out+"_scanvi_mod/training.pdf")
    
    data.obs["C_scANVI"] = lvae.predict(data)
    data.obsm["X_scANVI"] = lvae.get_latent_representation(data)
    
    sc.pp.neighbors(data, use_rep="X_scANVI")
    sc.tl.umap(data)
    data.obsm['X_scANVI_umap'] = data.obsm['X_umap']
    del data.obsm['X_umap']

    data.write_h5ad(args.h5ad_out)

    
if __name__ == "__main__":
    main()

