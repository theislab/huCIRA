import logging
import os

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _check_plot_deps() -> None:
    """Raise an informative error when optional plotting dependencies are missing."""
    missing = []
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        missing.append("matplotlib")
    try:
        import seaborn  # noqa: F401
    except ImportError:
        missing.append("seaborn")
    if missing:
        raise ImportError(
            f"Missing optional plotting dependencies: {', '.join(missing)}. "
            "Install them with: pip install 'hucira[plot]'"
        )


def _format_cytokine_names(x: str | list | np.ndarray | pd.Index) -> str | list[str]:
    if isinstance(x, (list, np.ndarray, pd.Index)):
        return [_format_cytokine_names(_x) for _x in x]
    text = x.get_text() if hasattr(x, "get_text") else x
    text = text.replace("beta", r"$\beta$")
    text = text.replace("alpha", r"$\alpha$")
    text = text.replace("gamma", r"$\gamma$")
    text = text.replace("lambda", r"$\lambda$")
    text = text.replace("omega", r"$\omega$")
    return text


def plot_significant_results(
    nes_plot_df: pd.DataFrame,
    annot_plot_df: pd.DataFrame,
    robust_results: pd.DataFrame | None = None,
    selected_celltypes: list[str] | None = None,
    selected_cytokines: list[str] | None = None,
    ax=None,
    fontsize: float = 6.0,
    save_fig: bool = False,
    fig_path: str = "",
    fig_width: float = 10.0,
    fig_height: float = 12.0,
) -> None:
    """Plot a heatmap of robust enrichment results for a single contrast.

    Parameters
    ----------
    nes_plot_df : pandas.DataFrame
        Pivot DataFrame of robust enrichment scores.
    annot_plot_df : pandas.DataFrame
        Annotation DataFrame of significance stars.
    robust_results : pandas.DataFrame or None
        Robust results DataFrame from :func:`~hucira.get_robust_significant_results`.
    selected_celltypes : list of str or None
        Subset of cell types to visualise.
    selected_cytokines : list of str or None
        Subset of cytokines to visualise.
    ax : matplotlib.axes.Axes or None
        An existing Axes to plot on.  If ``None``, a new figure is created.
    fontsize : float
        Font size for annotations and tick labels.
    save_fig : bool
        Whether to save the figure to disk.
    fig_path : str
        Directory for saved figures.
    fig_width : float
        Figure width in inches.
    fig_height : float
        Figure height in inches.
    """
    _check_plot_deps()

    import matplotlib.pyplot as plt
    import seaborn as sns

    if robust_results is not None and not robust_results.empty:
        contrast_name = robust_results.contrast.iloc[0]
    else:
        contrast_name = "Contrast1_vs_Contrast2"

    if isinstance(nes_plot_df, pd.DataFrame) and isinstance(annot_plot_df, pd.DataFrame):
        if selected_celltypes:
            nes_plot_df = nes_plot_df.T.loc[selected_celltypes].T
            annot_plot_df = annot_plot_df.T.loc[selected_celltypes].T
        if selected_cytokines:
            nes_plot_df = nes_plot_df.loc[selected_cytokines]
            annot_plot_df = annot_plot_df.loc[selected_cytokines]

        created_fig = ax is None
        if created_fig:
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        else:
            fig = ax.get_figure()
        sns.heatmap(
            nes_plot_df,
            square=True,
            annot=annot_plot_df,
            cmap="RdBu_r",
            center=0,
            annot_kws={"fontsize": fontsize, "family": "sans-serif"},
            fmt="",
            linewidths=0.5,
            linecolor="white",
            cbar=True,
            cbar_kws={"shrink": 0.5, "fraction": 0.04, "pad": 0.02},
            ax=ax,
        )
        ax.set_title(contrast_name, fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_facecolor("lightgray")
        ax.tick_params(axis="both", which="both", length=0)

        # Axis labels
        ax.set_xticks(0.5 + np.arange(nes_plot_df.shape[1]))
        ax.set_xticklabels(nes_plot_df.columns, fontsize=fontsize, rotation=90, ha="center")
        ax.set_yticks(0.5 + np.arange(nes_plot_df.shape[0]))
        ax.set_yticklabels(_format_cytokine_names(nes_plot_df.index), fontsize=fontsize, rotation=0, ha="right")

        if save_fig:
            fig.savefig(
                os.path.join(fig_path, f"heatmap_significant_results_{contrast_name}.svg"),
                bbox_inches="tight",
                pad_inches=0,
                dpi=500,
            )
        if created_fig:
            plt.show()
        return

    logger.warning("Nothing was plotted. Check input data!")
    return
