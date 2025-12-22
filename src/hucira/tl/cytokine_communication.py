import re
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData


warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


def _get_senders(
    adata: AnnData,
    cytokine_info: pd.DataFrame,
    cytokine: str = "IL-32-beta",
    show: bool = False,
    column_cell_type: str = "cell_type",
) -> pd.DataFrame:
    genes = np.unique(re.split(", ", cytokine_info.loc[cytokine_info.name == cytokine, "gene"].values[0]))
    mask = np.isin(genes, adata.var_names)

    if not mask.any():
        print(f"None of the cytokine producing genes ({genes}) were found in dataset for cytokine {cytokine}.")
        return None
    if not mask.all():
        print(f"The following cytokine producing genes were not found in the dataset and are excluded: {genes[~mask]}")
        genes = genes[mask]
    adata = adata[:, genes]
        
    # Ranks gene(s) of query sender cytokine across immune cell types.
    adata_out = sc.tl.rank_genes_groups(
        adata,
        groupby=column_cell_type,
        copy=True,
        use_raw=False,
        method="wilcoxon",
    )
    result = adata_out.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    results_mean, results_frac = [], []
    rank_genes_df = []
    for g in groups:
        df = pd.DataFrame(
            {
                "gene": result["names"][g],
                "logfoldchanges": result["logfoldchanges"][g],
                "pvals": result["pvals"][g],
                "pvals_adj": result["pvals_adj"][g],
                column_cell_type: g,
            }
        )
        rank_genes_df.append(df)
    rank_genes_df = pd.concat(rank_genes_df, axis=0)
    rank_genes_df.set_index(column_cell_type, inplace=True)
    grouped = rank_genes_df.groupby(column_cell_type)

    # Chooses minimum rank_genes_group() statistical parameters (considers limiting gene, if there are multiple per cytokine)
    grouped_rank_genes_df_all = []
    for celltype in grouped.groups.keys():
        grouped_celltype_df = grouped.get_group(celltype)

        # get gene with smallest log_fold_change (representing limiting gene), and retrieve stat. parameters
        limiting_gene_idx = np.argmin(grouped_celltype_df["logfoldchanges"].values)
        limiting_gene_vals = grouped_celltype_df.iloc[limiting_gene_idx][["logfoldchanges", "pvals", "pvals_adj"]]
        gene_concat = ", ".join(grouped_celltype_df["gene"])
        grouped_rank_genes_df = limiting_gene_vals.to_frame().T
        grouped_rank_genes_df["gene"] = gene_concat
        grouped_rank_genes_df.index = [celltype]
        grouped_rank_genes_df_all.append(grouped_rank_genes_df)

    grouped_rank_genes_df_all = pd.concat(grouped_rank_genes_df_all, axis=0)
    grouped_rank_genes_df_all = grouped_rank_genes_df_all.rename(
        columns={"logfoldchanges": "min_logfoldchanges", "pvals": "min_pvals", "pvals_adj": "min_pvals_adj"}
    )

    # Minimum of mean gene expression of sender cytokine genes:
    X_df = adata[:, genes].to_df()
    frac_df = X_df > 0
    X_df.loc[:, column_cell_type] = adata.obs.loc[:, column_cell_type].values
    frac_df.loc[:, column_cell_type] = adata.obs.loc[:, column_cell_type].values

    # take minimum average gene expression across all genes required for this sender
    results_mean = (
        X_df.groupby(column_cell_type, observed=False).mean().min(axis=1).to_frame().rename({0: "mean_X"}, axis=1)
    )
    # take minimum expression fraction across all genes required for this sender
    results_frac = (
        frac_df.groupby(column_cell_type, observed=False).mean().min(axis=1).to_frame().rename({0: "frac_X"}, axis=1)
    )

    # Final df with information about active sender cytokines.
    results = pd.concat([grouped_rank_genes_df_all, results_mean, results_frac], axis=1)
    results["mean_X>0"] = results["mean_X"].where(results["mean_X"] > 0, None)
    results.loc[:, "cytokine"] = cytokine
    return results
    

def _get_receivers(
    adata: AnnData, cytokine_info: pd.DataFrame, cytokine: str, column_cell_type: str = "cell_type"
) -> pd.DataFrame | None:
    # get receptor genes for this cytokine
    _receptor_genes = cytokine_info.loc[cytokine_info.name == cytokine, "receptor gene"]
    if _receptor_genes.isna().all():
        print(f"No receptor gene found in cytokine_info for cytokine: {cytokine}")
        return None
    assert len(_receptor_genes) == 1, _receptor_genes
    _receptor_genes = _receptor_genes.values[0]
    # there can be multiple receptors
    candidates = re.split("; ", _receptor_genes)
    results_mean, results_frac = [], []
    # each receptor may require the expression of multiple genes
    for candidate in candidates:
        # print(candidate)
        genes = np.array(re.split(", ", candidate))
        mask = np.isin(genes, adata.var_names)
        if not mask.any():
            print(f"None of the cytokine receptor genes ({genes}) were found in dataset for cytokine {cytokine}.")
            continue
        if not mask.all():
            print(
                f"The following cytokine receptor genes were not found in the dataset and are excluded: {genes[~mask]}"
            )
            genes = genes[mask]
        X_df = adata[:, genes].to_df()
        frac_df = X_df > 0
        X_df.loc[:, column_cell_type] = adata.obs.loc[:, column_cell_type].values
        frac_df.loc[:, column_cell_type] = adata.obs.loc[:, column_cell_type].values
        # take minimum average gene expression across all genes required for this receptor
        results_mean.append(X_df.groupby(column_cell_type, observed=False).mean().min(axis=1).to_frame())
        # take minimum expression fraction across all genes required for this receptor
        results_frac.append(frac_df.groupby(column_cell_type, observed=False).mean().min(axis=1).to_frame())
    if len(results_mean) == 0:
        return None

    results_mean = pd.concat(results_mean, axis=1).max(axis=1).to_frame().rename({0: "mean_X"}, axis=1)
    results_frac = pd.concat(results_frac, axis=1).max(axis=1).to_frame().rename({0: "frac_X"}, axis=1)
    results = pd.concat([results_mean, results_frac], axis=1)
    results.loc[:, "cytokine"] = cytokine
    return results


def get_one_senders_and_receivers(
    adata: AnnData,
    cytokine_info: pd.DataFrame,
    cytokine: str,
    celltype_colname: str = "cell_type",
    sender_pvalue_threshold: float = 0.1,
    receiver_mean_X_threshold: float = 0,
     sender_lfc_threshold: float = 0
) -> (pd.DataFrame, pd.DataFrame):
    """Generates cytokine producer and receiver statistics (senders and receivers of cell-cell communication) for one cytokine.

    Best for exploration purposes of a singular cytokine.

    Parameters
    ----------
    adata : AnnData
        Query adata object of analysis
    cytokine_info : pd.DataFrame
        External file containing info about receptor genes of each cytokine in format
        pd.DataFrame({"name": cytokine, "receptor gene": [gene1, gene2]})
    cytokine : str
        A cytokine, which ideally should be present in robust_results
        (the outcome of the robust enrichment analysis)
    celltype_colname : str, default "cell_type"
        Column name of where cell types are stored in adata

    Returns
    -------
    df_senders : pd.DataFrame
        Cytokine signal senders per cell type
    df_receivers : pd.DataFrame
        Cytokine signal receivers per cell type
    """
    df_senders = _get_senders(
        adata=adata, cytokine_info=cytokine_info, cytokine=cytokine, column_cell_type=celltype_colname
    )
    df_receivers = _get_receivers(
        adata=adata, cytokine_info=cytokine_info, cytokine=cytokine, column_cell_type=celltype_colname
    )
    if df_senders is not None:
        df_senders = df_senders.loc[
            (df_senders.min_pvals < sender_pvalue_threshold) & (df_senders.min_logfoldchanges > sender_lfc_threshold)
        ]
    if df_receivers is not None:
        df_receivers = df_receivers.loc[df_receivers.mean_X > receiver_mean_X_threshold]

    return df_senders, df_receivers


def get_all_senders_and_receivers(
    adata: AnnData,
    cytokine_info: pd.DataFrame,
    cytokine_list: list = None,
    celltype_colname: str = "cell_type",
    sender_pvalue_threshold: float = 0.1,
    receiver_mean_X_threshold: float = 0,
) -> (pd.DataFrame, pd.DataFrame):
    """Generates cytokine producer and receiver statistics (senders and receivers of cell-cell communication) for a list of cytokines.

    Best for visualization purposes (for plot_communication function).

    Parameters
    ----------
    adata : AnnData
        Query adata object of analysis
    cytokine_info : pd.DataFrame
        External file containing info about receptor genes of each cytokine in format
        pd.DataFrame({"name": cytokine, "receptor gene": [gene1, gene2]})
    cytokine_list : list, optional
        List of cytokines, which ideally should be present in robust_results
        (the outcome of the robust enrichment analysis). Default is None.
    celltype_colname : str, default "cell_type"
        Column name of where cell types are stored in adata

    Returns
    -------
    df_src : pd.DataFrame
        All cytokine signal senders
    df_tgt : pd.DataFrame
        All cytokine signal receivers
    """
    senders, receivers = [], []
    for cytokine in cytokine_list:
        df_senders, df_receivers = get_one_senders_and_receivers(
            adata=adata,
            cytokine_info=cytokine_info,
            cytokine=cytokine,
            celltype_colname=celltype_colname,
            sender_pvalue_threshold=0.1,
            receiver_mean_X_threshold=0,
        )

        if cytokine == "IL-32-beta":
            # no known receptor genes - create non-informative df_receivers manually.
            all_celltypes = sorted(adata.obs[celltype_colname].unique())
            df_receivers = pd.DataFrame.from_dict(
                dict(zip(all_celltypes, np.ones([len(all_celltypes), 2]) * np.inf, strict=True)),
                orient="index",
            ).rename({0: "mean_X", 1: "frac_X"}, axis=1)
            df_receivers.loc[:, "cytokine"] = cytokine

        if df_senders is not None and df_receivers is not None:
            df_senders = df_senders.assign(celltype=df_senders.index)
            df_receivers = df_receivers.assign(celltype=df_receivers.index)

            senders.append(df_senders)
            receivers.append(df_receivers)

    df_src = pd.concat(senders)
    df_tgt = pd.concat(receivers)

    return df_src, df_tgt
