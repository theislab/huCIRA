import logging
from typing import Literal

import gseapy as gp
import numpy as np
import pandas as pd
from anndata import AnnData

logger = logging.getLogger(__name__)


def _get_genesets(
    adata: AnnData,
    df: pd.DataFrame,
    celltype_signature: str,
    direction: Literal["upregulated", "downregulated", "both"] | None = None,
    threshold_pval: float | None = None,
    threshold_lfc: float | None = None,
) -> tuple[dict[str, list[str]], pd.DataFrame]:
    """Get shared gene sets between query adata and a reference gene-program DataFrame.

    Supports the Human Cytokine Dictionary, CIP signatures, or custom gene
    signatures for a chosen cell type.

    Parameters
    ----------
    adata : AnnData
        AnnData object with gene expression data.
    df : pandas.DataFrame
        Either the Human Cytokine Dictionary, CIP signatures, or a custom
        DataFrame containing columns ``["gene", "query_program", "celltype"]``.
    celltype_signature : str
        Cell-type name matching ``df.celltype``.
    direction : str or None
        One of ``"upregulated"``, ``"downregulated"``, or ``"both"``.
        Only relevant for the Human Cytokine Dictionary.
    threshold_pval : float or None
        Maximum adjusted p-value. Only relevant for the Human Cytokine
        Dictionary.
    threshold_lfc : float or None
        Minimum log-fold change. Only relevant for the Human Cytokine
        Dictionary.

    Returns
    -------
    gene_set_dict : dict
        Dictionary with cytokine/CIP as key and associated genes as values.
    gene_set_df : pandas.DataFrame
        DataFrame with gene overlap information between query data and gene
        program for the chosen cell type.
    """
    required_for_hcd = ["log_fc", "adj_p_value", "cytokine"]
    required_for_CIP = ["gene", "CIP", "celltype"]

    # Construct signature gene set if input is human cytokine dictionary
    if set(required_for_hcd).issubset(df.columns):
        logger.info("Computing gene sets of Human Cytokine Dictionary for %s.", celltype_signature)
        select = (df.adj_p_value <= threshold_pval) & (df.celltype == celltype_signature)
        if direction == "upregulated":
            select = select & (df.log_fc >= threshold_lfc)
        elif direction == "downregulated":
            select = select & (df.log_fc <= threshold_lfc)
        elif direction == "both":
            select = select & (df.log_fc.abs() >= threshold_lfc)
        else:
            raise ValueError(f"Invalid direction: {direction}.")
        df = df.loc[select]

        gene_set_dict = {}
        gene_set_df = pd.DataFrame()
        for cytokine_i, cytokine in enumerate(df.cytokine.unique()):
            gene_set = df.loc[df.cytokine == cytokine].gene.values
            gene_set_shared = np.intersect1d(gene_set, adata.var_names)
            gene_set_df.loc[cytokine_i, "cytokine"] = cytokine
            gene_set_df.loc[cytokine_i, "num_genes_signature"] = len(gene_set)
            gene_set_df.loc[cytokine_i, "num_shared_genes_signature"] = len(gene_set_shared)
            gene_set_df.loc[cytokine_i, "frac_shared_genes_signature"] = len(gene_set_shared) / len(gene_set)
            gene_set_dict[cytokine] = gene_set_shared

    # Construct signature gene set if input is CIP signatures
    elif set(required_for_CIP).issubset(df.columns):
        logger.info("Computing gene sets of Cytokine-induced gene programs for %s.", celltype_signature)
        select = df.celltype == celltype_signature
        df = df.loc[select]
        gene_set_dict = {}
        gene_set_df = pd.DataFrame()
        for CIP_i, CIP in enumerate(df.CIP.unique()):
            gene_set = df.loc[df.CIP == CIP].gene.values
            gene_set_shared = np.intersect1d(gene_set, adata.var_names)
            gene_set_df.loc[CIP_i, "CIP"] = CIP
            gene_set_df.loc[CIP_i, "num_genes_signature"] = len(gene_set)
            gene_set_df.loc[CIP_i, "num_shared_genes_signature"] = len(gene_set_shared)
            gene_set_df.loc[CIP_i, "frac_shared_genes_signature"] = len(gene_set_shared) / len(gene_set)
            gene_set_dict[CIP] = gene_set_shared

    # Construct signature gene set for custom gene programs
    elif "query_program" in df.columns:
        logger.info("Computing gene sets of user-defined gene programs for %s.", celltype_signature)
        select = df.celltype == celltype_signature
        df = df.loc[select]
        gene_set_dict = {}
        gene_set_df = pd.DataFrame()
        for query_program_i, query_program in enumerate(df.query_program.unique()):
            gene_set = df.loc[df.query_program == query_program].gene.values
            gene_set_shared = np.intersect1d(gene_set, adata.var_names)
            gene_set_df.loc[query_program_i, "query_program"] = query_program
            gene_set_df.loc[query_program_i, "num_genes_signature"] = len(gene_set)
            gene_set_df.loc[query_program_i, "num_shared_genes_signature"] = len(gene_set_shared)
            gene_set_df.loc[query_program_i, "frac_shared_genes_signature"] = len(gene_set_shared) / len(gene_set)
            gene_set_dict[query_program] = gene_set_shared

    else:
        raise ValueError(
            "invalid input for df parameter. You can use either the Human Cytokine Dictionary with load_human_cytokine_dict(), or our CIP signatures with load_CIP_signatures(). If you want to compute enrichment of custom gene sets, df must have columns: ['gene', 'query_program', 'celltype']."
        )
        return
    return gene_set_dict, gene_set_df


def _compute_mu_and_sigma(adata: AnnData, contrast_column: str, condition: str) -> dict:
    group = adata[adata.obs[contrast_column] == condition]
    num_cells = group.shape[0]
    X = group.X.toarray() if hasattr(group.X, "toarray") else group.X
    mu = np.mean(X, axis=0)
    sigma = np.std(X, axis=0, ddof=1)
    return {"mu": mu, "sigma": sigma, "num_cells": num_cells}


def _compute_s2n(
    adata: AnnData, contrast_column: str, condition_1: str, condition_2: str, precomputed_stats: dict | None = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute the signal-to-noise ratio (S2N) for each gene between two conditions.

    Parameters
    ----------
    adata : AnnData
        AnnData object with gene expression data.
    contrast_column : str
        Key in ``adata.obs`` indicating condition labels (e.g. ``"disease_state"``).
    condition_1 : str
        Name of the first condition (e.g. ``"flare"``).
    condition_2 : str
        Name of the second condition (e.g. ``"healthy"``).
    precomputed_stats : dict or None
        Pre-computed per-condition statistics from :func:`_compute_mu_and_sigma`.
        If provided, avoids recomputing mean and std.

    Returns
    -------
    stats : pandas.DataFrame
        S2N values indexed by gene names.
    num_cells : pandas.DataFrame
        Number of cells per condition.
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
        logger.debug("Using precomputed stats")
        num_cells_1 = precomputed_stats[condition_1]["num_cells"]
        num_cells_2 = precomputed_stats[condition_2]["num_cells"]
        mu1 = precomputed_stats[condition_1]["mu"]
        mu2 = precomputed_stats[condition_2]["mu"]
        sigma1 = precomputed_stats[condition_1]["sigma"]
        sigma2 = precomputed_stats[condition_2]["sigma"]

    # Compute S2N
    s2n = (mu1 - mu2) / (sigma1 + sigma2 + 1e-8)  # epsilon to avoid division by zero

    num_cells = pd.DataFrame(
        index=[f"{condition_1}_vs_{condition_2}"],
        columns=["num_cells_1", "num_cells_2"],
        data=[[num_cells_1, num_cells_2]],
    )
    stats = pd.DataFrame(s2n, index=adata.var_names, columns=[f"{condition_1}_vs_{condition_2}"])

    return stats, num_cells


def _compute_ranking_statistic(
    adata: AnnData, contrast_column: str, contrasts_combo: list[tuple[str, str]]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rnk_stats, num_cells = [], []
    precomputed_stats = {}

    conditions = []
    for condition in contrasts_combo:
        conditions.extend([condition[0], condition[1]])
    conditions = np.unique(conditions)

    for condition in conditions:
        precomputed_stats[condition] = _compute_mu_and_sigma(
            adata, contrast_column=contrast_column, condition=condition
        )

    for condition in contrasts_combo:
        _rnk_stats, _num_cells = _compute_s2n(
            adata,
            contrast_column=contrast_column,
            condition_1=condition[0],
            condition_2=condition[1],
            precomputed_stats=precomputed_stats,
        )
        rnk_stats.append(_rnk_stats)
        num_cells.append(_num_cells)
    return pd.concat(rnk_stats, axis=1), pd.concat(num_cells, axis=0)


def run_one_enrichment_test(
    adata: AnnData,
    df: pd.DataFrame,
    celltype_combo: tuple[str, str] = ("B cell", "B_cell"),
    celltype_column: str = "cell_type",
    contrasts_combo: tuple[str, str] | list[tuple[str, str]] = None,
    contrast_column: str = "disease_state",
    direction: Literal["upregulated", "downregulated", "both"] = "upregulated",
    # Filtering parameters for gene set construction
    threshold_lfc: float = 1.0,
    threshold_expression: float = 0.0,
    threshold_pval: float = 0.01,
    # GSEA parameters
    min_size: int = 10,
    max_size: int = 1000,
    permutation_num: int = 1000,
    weight: float = 1.0,
    seed: int = 2025,
    verbose: bool = False,
    threads: int = 6,
) -> pd.DataFrame | None:
    """Compute cytokine enrichment activity in one cell type using GSEA scoring.

    1. "Looks up" query cell type in human cytokine dictionary and retrieves
       associated up-/downregulated genes per cytokine as reference.
    2. Creates ranking of query data genes contrasting condition 1 vs
       condition 2.
    3. Computes enrichment of each cytokine by matching their associated gene
       set in the ranked list.

    Parameters
    ----------
    adata : ~anndata.AnnData
        The query AnnData object.
    df : pandas.DataFrame
        Human Cytokine Dictionary, CIP signatures, or custom gene-program
        DataFrame.
    celltype_combo : tuple of (str, str)
        Cell-type name in the query *adata* (first) and the corresponding
        cell-type name in *df* (second).
    celltype_column : str
        Column name in ``adata.obs`` that stores the cell types.
    contrasts_combo : tuple of (str, str) or list of tuple
        Pair(s) of biological conditions to compare.
    contrast_column : str
        Column name in ``adata.obs`` that stores condition labels.
    direction : str
        One of ``"upregulated"``, ``"downregulated"``, or ``"both"``.
    threshold_lfc : float
        Minimum log-fold change used to build the gene set.
    threshold_expression : float
        Minimum mean gene expression across all cells.
    threshold_pval : float
        Maximum adjusted p-value used to build the gene set.
    min_size : int
        Minimum gene-set size passed to GSEA.
    max_size : int
        Maximum gene-set size passed to GSEA.
    permutation_num : int
        Number of permutations for GSEA.
    weight : float
        GSEA weighting exponent.
    seed : int
        Random seed for reproducibility.
    verbose : bool
        Whether to print progress messages.
    threads : int
        Number of threads for GSEA.

    Returns
    -------
    pandas.DataFrame
        DataFrame with all computed enrichment scores and statistical
        parameters.  Not filtered by significance or robustness yet.
    """
    if not isinstance(contrasts_combo, list):
        assert isinstance(contrasts_combo, tuple)
        contrasts_combo = [contrasts_combo]

    celltype_adata = celltype_combo[0]
    celltype_signature = celltype_combo[1]

    # allows potential loop of celltype combos to continue
    if celltype_adata not in adata.obs[celltype_column].unique():
        logger.warning(
            "'%s' is not present in celltype_column (%s) of query adata. Skipping enrichment test of this celltype.",
            celltype_adata,
            celltype_column,
        )
        return None

    # filter for cell type
    logger.debug("Filter for cell type:")
    adata = adata[adata.obs[celltype_column] == celltype_adata]
    logger.debug("Filter for cell type: done.")

    # filter based on gene expression
    logger.debug("Filter for gene expression:")
    adata = adata[:, adata.X.mean(axis=0) >= threshold_expression]
    logger.debug("Filter for gene expression: done.")

    # get genesets
    logger.debug("Get gene sets:")
    gene_set_dict, gene_set_df = _get_genesets(
        adata=adata,
        df=df,
        celltype_signature=celltype_signature,
        direction=direction,
        threshold_pval=threshold_pval,
        threshold_lfc=threshold_lfc,
    )

    logger.debug("Get gene sets: done.")

    # compute ranking stat
    logger.debug("Get ranking stats:")
    rnk_stats, num_cells_per_condition = _compute_ranking_statistic(
        adata, contrast_column=contrast_column, contrasts_combo=contrasts_combo
    )
    logger.debug("Get ranking stats: done.")
    results = []

    for contrast_name in rnk_stats.columns:
        logger.info("Running enrichment for contrast: %s", contrast_name)
        # format stat so that it can be processed be gseapy
        rnk = (
            rnk_stats.loc[:, contrast_name]
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .sort_values(ascending=False)
            .to_frame()
            .rename({contrast_name: 1}, axis=1)
            .rename_axis("0")
        )

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
        logger.debug("%s: done.", contrast_name)

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

    required_for_hcd = ["log_fc", "adj_p_value", "cytokine"]
    if set(required_for_hcd).issubset(df.columns):
        results.rename({"Term": "cytokine"}, inplace=True, axis=1)
        results = pd.merge(results, gene_set_df, on="cytokine")
    elif "CIP" in df.columns:
        results.rename({"Term": "CIP"}, inplace=True, axis=1)
        results = pd.merge(results, gene_set_df, on="CIP")
        results.direction = "upregulated"
    elif "query_program" in df.columns:
        results.rename({"Term": "query_program"}, inplace=True, axis=1)
        results = pd.merge(results, gene_set_df, on="query_program")
        results.direction = "custom input"

    return results


def run_all_enrichment_test(
    adata: AnnData,
    df: pd.DataFrame,
    celltype_combos: list[tuple[str, str]] = None,
    celltype_column: str = "cell_type",
    contrasts_combo: tuple[str, str] | list[tuple[str, str]] = None,
    contrast_column: str = "disease_state",
    direction: Literal["upregulated", "downregulated", "both"] = "upregulated",
    # Filtering parameters for gene set construction
    threshold_lfc: float | list[float] = 1.0,
    threshold_expression: float | list[float] = 0.0,
    threshold_pval: float = 0.01,
    # GSEA parameters
    min_size: int = 10,
    max_size: int = 1000,
    permutation_num: int = 1000,
    weight: float = 1.0,
    seed: int = 2025,
    verbose: bool = False,
    threads: int = 6,
) -> pd.DataFrame:
    """Compute cytokine enrichment across multiple threshold values for robustness.

    Wrapper around :func:`run_one_enrichment_test` that loops through several
    threshold values to obtain more robust results.

    Parameters
    ----------
    adata : ~anndata.AnnData
        The query AnnData object.
    df : pandas.DataFrame
        Human Cytokine Dictionary, CIP signatures, or custom gene-program
        DataFrame.
    celltype_combos : list of tuple
        List of ``(query_celltype, df_celltype)`` pairs.
    celltype_column : str
        Column name in ``adata.obs`` that stores the cell types.
    contrasts_combo : tuple of (str, str) or list of tuple
        Pair(s) of biological conditions to compare.
    contrast_column : str
        Column name in ``adata.obs`` that stores condition labels.
    direction : str
        One of ``"upregulated"``, ``"downregulated"``, or ``"both"``.
    threshold_lfc : float or list of float
        Log-fold change threshold(s) for building the gene set.
    threshold_expression : float or list of float
        Mean gene-expression threshold(s).
    threshold_pval : float
        Maximum adjusted p-value for the gene set.
    min_size : int
        Minimum gene-set size passed to GSEA.
    max_size : int
        Maximum gene-set size passed to GSEA.
    permutation_num : int
        Number of permutations for GSEA.
    weight : float
        GSEA weighting exponent.
    seed : int
        Random seed for reproducibility.
    verbose : bool
        Whether to print progress messages.
    threads : int
        Number of threads for GSEA.

    Returns
    -------
    pandas.DataFrame
        DataFrame with all computed enrichment scores and statistical
        parameters from multiple thresholds.
    """
    if isinstance(threshold_lfc, float):
        threshold_lfc = [threshold_lfc]
    if isinstance(threshold_expression, float):
        threshold_expression = [threshold_expression]

    all_enrichment_results = []
    for _celltype_combo_k, celltype_combo in enumerate(celltype_combos):
        for lfc in threshold_lfc:
            for expr in threshold_expression:
                results = run_one_enrichment_test(
                    adata=adata,
                    df=df,
                    celltype_combo=celltype_combo,
                    celltype_column=celltype_column,
                    contrasts_combo=contrasts_combo,
                    contrast_column=contrast_column,
                    direction=direction,
                    # Robustness parameters
                    threshold_pval=threshold_pval,
                    threshold_lfc=lfc,
                    threshold_expression=expr,
                    # GSEA parameters
                    min_size=min_size,
                    max_size=max_size,
                    permutation_num=permutation_num,
                    weight=weight,
                    seed=seed,
                    verbose=verbose,
                    threads=threads,
                )

                all_enrichment_results.append(results)

    all_enrichment_results = pd.concat(all_enrichment_results, axis=0)

    return all_enrichment_results
