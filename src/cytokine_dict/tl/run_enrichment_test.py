from anndata import AnnData
import numpy as np
import pandas as pd
# import scanpy as sc
import seaborn as sns
import gseapy as gp
from typing import Tuple, Literal


def vprint(msg, verbose):
    if verbose:
        print(msg)


def get_genesets(
    adata: AnnData,
    df_hcd_all: pd.DataFrame,
    celltype_signature: str,
    direction: Literal["upregulated", "downregulated", "both"],
    threshold_pval: float,
    threshold_lfc: float, 
) -> Tuple[dict[str, list[str]], pd.DataFrame]:

    # construct signature gene set
    select = ((df_hcd_all.adj_p_value <= threshold_pval) & (df_hcd_all.celltype == celltype_signature))
    if direction == "upregulated":
        select = select & (df_hcd_all.log_fc >= threshold_lfc)
    elif direction == "downregulated":
        select = select & (df_hcd_all.log_fc <= threshold_lfc)
    elif direction == "both":
        select = select & (df_hcd_all.log_fc.abs() >= threshold_lfc)         
    else:
        raise ValueError(f"Invalid direction: {direction}.")
    
    df_hcd_ct = df_hcd_all.loc[select]
    
    gene_set_dict = {}
    gene_set_df = pd.DataFrame()
    for cytokine_i, cytokine in enumerate(df_hcd_ct.cytokine.unique()):
        gene_set = df_hcd_ct.loc[df_hcd_ct.cytokine == cytokine].gene.values
        gene_set_shared = np.intersect1d(gene_set, adata.var_names)
        gene_set_df.loc[cytokine_i, "cytokine"] = cytokine
        gene_set_df.loc[cytokine_i, "num_genes_signature"] = len(gene_set)
        gene_set_df.loc[cytokine_i, "num_shared_genes_signature"] = len(gene_set_shared)
        gene_set_df.loc[cytokine_i, "frac_shared_genes_signature"] = len(gene_set_shared) / len(gene_set)
        gene_set_dict[cytokine] = gene_set_shared        

    return gene_set_dict, gene_set_df


def compute_mu_and_sigma(adata: AnnData, contrast_column: str, condition: str) -> pd.DataFrame:

    group = adata[adata.obs[contrast_column] == condition]
    num_cells = group.shape[0]
    X = group.X.toarray() if hasattr(group.X, "toarray") else group.X
    mu = np.mean(X, axis=0)
    sigma = np.std(X, axis=0, ddof=1)
    return {"mu": mu, "sigma": sigma, "num_cells": num_cells}


def compute_s2n(
    adata: AnnData, 
    contrast_column: str, 
    condition_1: str, 
    condition_2: str,
    precomputed_stats: dict | None = None
    
) -> (pd.DataFrame, pd.DataFrame):
    """
    Compute the signal-to-noise ratio (S2N) for each gene between two conditions in an AnnData object.

    Parameters:
    - adata: AnnData object with gene expression data.
    - contrast_column: Key in `adata.obs` indicating the condition labels (e.g. "disease_state").
    - condition_1: Name of the first condition (e.g., "flare").
    - condition_2: Name of the second condition (e.g., "healthy").
    
    Returns:
    - s2n_scores: pandas Series of S2N values indexed by gene names.
    """

    if precomputed_stats is None:
    
        # Select cells for each condition
        group1 = adata[adata.obs[contrast_column] == condition_1]
        group2 = adata[adata.obs[contrast_column] == condition_2]
        
        # number of cells per condition
        num_cells_1 = group1.shape[0]
        num_cells_2 = group2.shape[0]
        
        # Get expression matrices
        X1 = group1.X.toarray() if hasattr(group1.X, "toarray") else group1.X
        X2 = group2.X.toarray() if hasattr(group2.X, "toarray") else group2.X
    
        # Compute mean and std per gene
        mu1 = np.mean(X1, axis=0)
        mu2 = np.mean(X2, axis=0)
        sigma1 = np.std(X1, axis=0, ddof=1)
        sigma2 = np.std(X2, axis=0, ddof=1)

    else:
        vprint("Using precomputed stats", True)
        num_cells_1 = precomputed_stats[condition_1]["num_cells"]
        num_cells_2 = precomputed_stats[condition_2]["num_cells"]
        mu1 = precomputed_stats[condition_1]["mu"]
        mu2 = precomputed_stats[condition_2]["mu"]
        sigma1 = precomputed_stats[condition_1]["sigma"]
        sigma2 = precomputed_stats[condition_2]["sigma"]

    # Compute S2N
    s2n = (mu1 - mu2) / (sigma1 + sigma2 + 1e-8)  # epsilon to avoid division by zero

    num_cells = pd.DataFrame(index=[f"{condition_1}_vs_{condition_2}"], columns=["num_cells_1", "num_cells_2"], data=[[num_cells_1, num_cells_2]])
    stats = pd.DataFrame(s2n, index=adata.var_names, columns=[f"{condition_1}_vs_{condition_2}"])
    
    return stats, num_cells

  
def compute_ranking_statistic(adata: AnnData, contrast_column: str, contrasts: list[Tuple[str, str]]) -> (pd.DataFrame, pd.DataFrame):
    rnk_stats, num_cells = [], []
    precomputed_stats = {}

    conditions = []
    for condition in contrasts:
        conditions.extend([condition[0], condition[1]])
    conditions = np.unique(conditions)
    
    for condition in conditions:
        precomputed_stats[condition] = compute_mu_and_sigma(
            adata, 
            contrast_column=contrast_column, 
            condition=condition
        )

    for condition in contrasts:
        _rnk_stats, _num_cells = compute_s2n(
            adata,
            contrast_column=contrast_column,
            condition_1=condition[0],
            condition_2=condition[1],
            precomputed_stats=precomputed_stats,
        )
        rnk_stats.append(_rnk_stats)
        num_cells.append(_num_cells)
    return pd.concat(rnk_stats, axis=1), pd.concat(num_cells, axis=0)

    
def run_enrichment_test(
    adata: AnnData,
    df_hcd_all: pd.DataFrame,
    celltype: Tuple[str, str] = ("B cell", "B"),
    direction: str = "upregulated",
    threshold_pval: float = 0.01,
    threshold_lfc: float = 1.,
    threshold_expression: float = 0.0,
    contrast_column: str = "disease_state",
    celltype_column: str = "disease_state",
    contrasts: Tuple[str, str] | list[Tuple[str, str]] = None,
    min_size: int = 10,
    max_size: int = 1000,
    permutation_num: int = 1000,
    weight: float = 1.,
    seed: int = 2025,
    verbose: bool = True,
    threads: int = 6,
) -> pd.DataFrame:

    if not isinstance(contrasts, list):
        assert isinstance(contrasts, Tuple)
        contrasts = [contrasts]
        
    celltype_adata = celltype[0]
    celltype_signature = celltype[1]

    # allows potential loop of celltype combos to continue
    if celltype_adata not in adata.obs[celltype_column].unique():
        print(f"'{celltype_adata}' is not present in celltype_column ({celltype_column}) of query adata. Skipping enrichment test of this celltype.\n") 
        return None

    # filter for cell type
    vprint("Filter for cell type:", verbose)
    adata = adata[adata.obs[celltype_column] == celltype_adata]
    vprint("Filter for cell type: done.", verbose)
    
    # filter based on gene expression
    vprint("Filter for gene expression:", verbose)
    adata = adata[:, adata.X.mean(axis=0) >= threshold_expression]
    vprint("Filter for gene expression: done.", verbose)

    # get genesets
    vprint("Get gene sets:", verbose)
    gene_set_dict, gene_set_df = get_genesets(
        adata=adata,
        df_hcd_all=df_hcd_all,
        celltype_signature=celltype_signature,
        direction=direction,
        threshold_pval=threshold_pval,
        threshold_lfc=threshold_lfc,
    )

    vprint("Get gene sets: done.", verbose)
    
    # compute ranking stat
    vprint("Get ranking stats:", verbose)
    rnk_stats, num_cells_per_condition = compute_ranking_statistic(adata, contrast_column=contrast_column, contrasts=contrasts)
    vprint("Get ranking stats: done.", verbose)
    results = []
    
    for contrast_name in rnk_stats.columns:
        print(contrast_name)
        # format stat so that it can be processed be gseapy
        rnk = rnk_stats.loc[:, contrast_name]\
            .replace([np.inf, -np.inf], np.nan)\
            .dropna()\
            .sort_values(ascending=False)\
            .to_frame()\
            .rename({contrast_name: 1}, axis=1)\
            .rename_axis("0")
        
        # run enrichment
        gp_res = gp.prerank(
            rnk=rnk,
            gene_sets=gene_set_dict,
            min_size=min_size,
            max_size=max_size,
            permutation_num=permutation_num,
            weight=weight,
            outdir=None,
            seed=seed,
            verbose=verbose,
            threads=threads,
        )
        _res = gp_res.res2d
        _res.loc[:, "contrast"] = contrast_name
        _res.loc[:, "num_cells_1"] = num_cells_per_condition.loc[contrast_name, "num_cells_1"]
        _res.loc[:, "num_cells_2"] = num_cells_per_condition.loc[contrast_name, "num_cells_2"]
        _res.loc[:, "percent_duplicate_ranking_stats"] = (rnk.duplicated(keep="first").sum() / rnk.shape[0]) * 100
        results.append(_res)
        vprint(f"{contrast_name}: done.", verbose)

    # combine results and save hyperparams
    results = pd.concat(results, axis=0, ignore_index=True)
    results.loc[:, "celltype_adata"] = celltype_adata
    results.loc[:, "celltype_signature"] = celltype_signature
    results.loc[:, "celltype_combo"] = f"{celltype_adata} ({celltype_signature})"
    results.loc[:, "direction"] = direction
    results.loc[:, "threshold_pval"] = threshold_pval
    results.loc[:, "threshold_lfc"] = threshold_lfc
    results.loc[:, "threshold_expression"] = threshold_expression
    results.loc[:, "min_size"] = min_size
    results.loc[:, "max_size"] = max_size
    results.loc[:, "permutation_num"] = permutation_num
    results.loc[:, "weight"] = weight
    results.loc[:, "seed"] = seed
    results.loc[:, "threads"] = threads

    results.rename({"Term": "cytokine"}, inplace=True, axis=1)
    results = pd.merge(results, gene_set_df, on="cytokine")
    
    return results
