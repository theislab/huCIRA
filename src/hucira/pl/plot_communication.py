import re
import warnings

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from bokeh.palettes import all_palettes
from pycirclize import Circos

warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

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
        cytokine_colors = cytokine_colors[:len(unique_cytokines)] # in case there are less than 3 unique cytokines
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
