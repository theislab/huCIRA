import re
import warnings

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from bokeh.palettes import all_palettes
from pycirclize import Circos

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

    # Ranks genes of sender cytokine across immune cell types.
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

        # get gene with smallest abs(log_fold_change), representing limiting gene
        limiting_gene_idx = grouped_celltype_df["logfoldchanges"].idxmin()
        limiting_gene_vals = grouped_celltype_df.loc[limiting_gene_idx, ["logfoldchanges", "pvals", "pvals_adj"]]

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

    # Final df with informationa about active sender cytokines.
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
            (df_senders.min_pvals < sender_pvalue_threshold) & (df_senders.min_logfoldchanges > 0)
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


def plot_communication(
    df_src: pd.DataFrame,
    df_tgt: pd.DataFrame,
    frac_expressing_cells_sender: float | None = 0.05,
    frac_expressing_cells_receiver: float | None = 0.05,
    mean_cytokine_gene_expression_sender: float | None = None,
    mean_cytokine_gene_expression_receiver: float | None = None,
    df_enrichment: pd.DataFrame | None = None,
    all_celltypes: list | None = None,
    cytokine2color: dict | None = None,
    celltype2color: dict | None = None,
    figsize: tuple[float, float] = (5, 5),
    show_legend: bool = True,
    save_path: str | None = None,
    lw: float = 1.0,
    fontsize: int = 6,
    loc: str = "upper left",
    bbox_to_anchor: tuple[float, float] = (1, 1),
):
    """Generates a Circos plot to visualize cell-cell communication based on cytokine producers and receivers.

    The function filters the input dataframes based on thresholds for fraction of expressing cells
    and mean cytokine gene expression, then creates a circular layout with cell type partitions
    and draws directed links representing cytokine communication between producers and receivers.

    Parameters
    ----------
    df_src : pd.DataFrame
        DataFrame containing producer cell type and cytokine expression statistics,
        typically from `_get_expression_stats`. Must have 'celltype', 'cytokine',
        'mean_cytokine_gene_expression', and 'frac_expressing_cells' columns.
    df_tgt : pd.DataFrame
        DataFrame containing receiver cell type and cytokine expression statistics,
        typically from `_get_expression_stats`. Must have 'celltype', 'cytokine',
        'mean_cytokine_gene_expression', and 'frac_expressing_cells' columns.
    frac_expressing_cells_sender : float | None, default 0.05
        Minimum fraction of cells expressing a cytokine gene for a producer cell type.
        If None, no filtering is applied.
    frac_expressing_cells_receiver : float | None, default 0.05
        Minimum fraction of cells expressing a cytokine gene for a receiver cell type.
        If None, no filtering is applied.
    mean_cytokine_gene_expression_sender : float | None, default None
        Minimum mean expression of a cytokine gene for a producer cell type. If None, no filtering is applied.
    mean_cytokine_gene_expression_receiver : float | None, default None
        Minimum mean expression of a cytokine gene for a receiver cell type. If None, no filtering is applied.
    df_enrichment : pd.DataFrame | None, optional
        Optional dataframe with enrichment information. Default is None.
    all_celltypes : list | None, optional
        List of all cell types. If None, inferred from df_src and df_tgt.
    cytokine2color : dict | None, optional
        Optional mapping from cytokine names to colors.
    celltype2color : dict | None, optional
        Optional mapping from cell type names to colors.
    figsize : tuple[float, float], default (5, 5)
        Figure size for the plot.
    show_legend : bool, default True
        Whether to show the legend.
    save_path : str | None, optional
        Path to save the figure. If None, figure is not saved.
    lw : float, default 1.0
        Line width for links.
    fontsize : int, default 6
        Font size for labels.
    loc : str, default "upper left"
        Legend location.
    bbox_to_anchor : tuple[float, float], default (1, 1)
        Bounding box anchor for the legend.

    """
    if frac_expressing_cells_sender is not None:
        df_src = df_src.loc[df_src.frac_X > frac_expressing_cells_sender]
    if frac_expressing_cells_receiver is not None:
        df_tgt = df_tgt.loc[df_tgt.frac_X > frac_expressing_cells_receiver]
    if mean_cytokine_gene_expression_sender is not None:
        df_src = df_src.loc[df_src.mean_X > mean_cytokine_gene_expression_sender]
    if frac_expressing_cells_receiver is not None:
        df_tgt = df_tgt.loc[df_tgt.mean_X > mean_cytokine_gene_expression_receiver]

    if all_celltypes is None:
        all_celltypes = sorted(np.union1d(df_src.celltype.unique(), df_tgt.celltype.unique()))
    # celltype_colors = all_palettes["Set3"][len(all_celltypes)]
    if celltype2color is None:
        n = len(all_celltypes)

        # Get first 20 colors from Category20
        palette_20 = all_palettes["Category20"][20]
        # Get 20 colors from Category20b
        palette_20b = all_palettes["Category20b"][20]

        # Combine palettes
        combined_palette = palette_20 + palette_20b

        if n > 40:
            raise ValueError(f"Too many cell types ({n}) for available palettes (max 40).")

        # Assign colors to cell types
        celltype_colors = combined_palette[:n]
        celltype2color = dict(zip(all_celltypes, celltype_colors, strict=True))

    all_cytokines = np.union1d(df_src.cytokine.unique(), df_tgt.cytokine.unique())
    cytokine2idx = {cytokine: k for k, cytokine in enumerate(all_cytokines)}
    # cytokine_colors = all_palettes["Category20"][len(all_cytokines)]
    # cytokine2color = dict(zip(all_cytokines, cytokine_colors, strict=True))

    unique_cytokines = df_src.cytokine.unique()
    if df_enrichment is not None:
        significant_cytokines = df_enrichment.cytokine.unique()
        unique_cytokines = np.intersect1d(unique_cytokines, significant_cytokines)

    if cytokine2color is None:
        cytokine_colors = all_palettes["Colorblind"][max(3, len(unique_cytokines))]
        # cytokine_colors = all_palettes["Set3"][max(3, len(unique_cytokines))]
        cytokine2color = dict(zip(unique_cytokines, cytokine_colors, strict=True))

    # draw outer circle / cell type partitions
    sectors = dict(zip(all_celltypes, (2 * len(all_cytokines) + 3) * np.ones(len(all_celltypes)), strict=True))

    circos = Circos(sectors, space=3)
    for sector in circos.sectors:
        start, stop = sector.deg_lim
        center = (start + stop) / 2
        track = sector.add_track((92, 100))

        if 160 >= center >= 20:
            ha = "left"
        elif 340 >= center >= 200:
            ha = "right"
        else:
            ha = "center"

        if center < 90 or center > 270:
            va = "bottom"
        else:
            va = "top"

        track.axis(facecolor=celltype2color[sector.name])
        # track.text(shorten_cell_type_names(sector.name), color="black", size=6, r=110, rotation="horizontal", adjust_rotation=False, family="sans-serif", ha=ha)
        track.text(
            sector.name,
            color="black",
            size=fontsize,
            r=110,
            rotation="horizontal",
            adjust_rotation=False,
            family="sans-serif",
            ha=ha,
            va=va,
        )

    # draw links
    legend_cytokine2color = {}
    for _row_idx, row in df_src.iterrows():
        src_celltype = row.celltype
        cytokine_idx = cytokine2idx[row.cytokine]
        tgt_celltypes = df_tgt.loc[df_tgt.cytokine == row.cytokine].celltype.unique()

        for tgt_celltype in tgt_celltypes:
            is_enriched = True  # default --> plot if enriched or whenever no enrichment info is provided

            if df_enrichment is not None:
                df_enrichment.loc[:, "celltype"] = df_enrichment.celltype_combo.apply(lambda x: x.split(" (")[0])
                select = (df_enrichment.celltype == tgt_celltype) & (df_enrichment.cytokine == row.cytokine)
                is_enriched = df_enrichment.loc[select].shape[0] > 0

            if is_enriched:
                linestyle = None
                _score = df_tgt.loc[
                    (df_tgt.cytokine == row.cytokine) & (df_tgt.celltype == tgt_celltype), "mean_X"
                ].values
                assert len(_score) == 1
                if not np.isfinite(_score[0]):
                    linestyle = "--"

                circos.link_line(
                    (src_celltype, 1 + cytokine_idx),  # src node
                    (tgt_celltype, 2 + len(all_cytokines) + cytokine_idx),  # tgt node
                    direction=1,
                    color=cytokine2color[row.cytokine],
                    # color=celltype2color[src_celltype],
                    lw=lw,
                    arrow_height=8.0,
                    arrow_width=8.0,
                    linestyle=linestyle,
                )
                if row.cytokine not in legend_cytokine2color.keys():
                    legend_cytokine2color[row.cytokine] = cytokine2color[row.cytokine]

    circos.plotfig(figsize=figsize)
    plt.gca()

    legend_handles = []
    legend_labels = []
    for cytokine, color in legend_cytokine2color.items():
        legend_handles.append(mlines.Line2D([], [], color=color, lw=1.5))
        legend_labels.append(cytokine)
    if show_legend:
        plt.legend(
            handles=legend_handles,
            labels=legend_labels,
            title="Cytokines",
            loc=loc,
            bbox_to_anchor=bbox_to_anchor,
            prop={"family": "sans-serif", "size": 6},
            title_fontsize=6,
        )
    plt.tight_layout()
    if save_path:
        plt.savefig(
            save_path,
            bbox_inches="tight",
            pad_inches=0,
            transparent=True,
            dpi=400,
        )
    plt.show()

    return legend_handles, legend_labels
