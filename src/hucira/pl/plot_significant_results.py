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
    results_pivot: pd.DataFrame,
    df_annot: pd.DataFrame,
    robust_results_dict: dict[str, pd.DataFrame] | None = None,
    selected_celltypes: list[str] | None = None,
    selected_cytokines: list[str] | None = None,
    fontsize: float = 6.0,
    save_fig: bool = False,
    fig_path: str = "",
    fig_width: float = 10.0,
    fig_height: float = 12.0,
) -> None:
    """Plot a heatmap of robust enrichment results.

    Plots either the robust results from a dict of contrasts or individually
    per contrast.

    Parameters
    ----------
    results_pivot : pandas.DataFrame
        Pivot DataFrame of robust enrichment scores for one contrast.
    df_annot : pandas.DataFrame
        Annotation DataFrame of significance stars for one contrast.
    robust_results_dict : dict or None
        Robust enrichment score dictionary from
        :func:`~hucira.get_robust_significant_results`.  If provided, takes
        precedence over *results_pivot* and *df_annot*.
    selected_celltypes : list of str or None
        Subset of cell types to visualise.
    selected_cytokines : list of str or None
        Subset of cytokines to visualise.
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

    # Case 1: robust_results_dict is provided. This precedes the other arguments.
    if robust_results_dict is not None and len(robust_results_dict) > 0:
        n = len(robust_results_dict)
        fig, axes = plt.subplots(1, n, squeeze=False)

        for i, (contrast, (pivot, annot, _)) in enumerate(robust_results_dict.items()):
            ax = axes[0, i]

            # Apply filtering if requested
            if selected_celltypes:
                pivot = pivot.T.loc[selected_celltypes].T
                annot = annot.T.loc[selected_celltypes].T
            if selected_cytokines:
                pivot = pivot.loc[selected_cytokines]
                annot = annot.loc[selected_cytokines]

            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            sns.heatmap(
                pivot,
                square=True,
                annot=annot,
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

            ax.set_title(contrast, fontsize=10)
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.set_facecolor("lightgray")
            ax.tick_params(axis="both", which="both", length=0)

            # Axis labels
            ax.set_xticks(0.5 + np.arange(pivot.shape[1]))
            ax.set_xticklabels(pivot.columns, fontsize=fontsize, rotation=90, ha="center")
            ax.set_yticks(0.5 + np.arange(pivot.shape[0]))
            ax.set_yticklabels(_format_cytokine_names(pivot.index), fontsize=fontsize, rotation=0, ha="right")

        if save_fig:
            plt.savefig(
                os.path.join(fig_path, "all_contrasts_significant_results.svg"),
                bbox_inches="tight",
                pad_inches=0,
                dpi=500,
            )
        plt.tight_layout()
        plt.show()
        return

    # Case 2: single robust_result is provided, only the one chosen contrast comparison is plotted.
    if isinstance(results_pivot, pd.DataFrame) and isinstance(df_annot, pd.DataFrame):
        if selected_celltypes:
            results_pivot = results_pivot.T.loc[selected_celltypes].T
            df_annot = df_annot.T.loc[selected_celltypes].T
        if selected_cytokines:
            results_pivot = results_pivot.loc[selected_cytokines]
            df_annot = df_annot.loc[selected_cytokines]

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        sns.heatmap(
            results_pivot,
            square=True,
            annot=df_annot,
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
        ax.set_title("Contrast1_vs_Contrast2", fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_facecolor("lightgray")
        ax.tick_params(axis="both", which="both", length=0)

        # Axis labels
        ax.set_xticks(0.5 + np.arange(results_pivot.shape[1]))
        ax.set_xticklabels(results_pivot.columns, fontsize=fontsize, rotation=90, ha="center")
        ax.set_yticks(0.5 + np.arange(results_pivot.shape[0]))
        ax.set_yticklabels(_format_cytokine_names(results_pivot.index), fontsize=fontsize, rotation=0, ha="right")

        plt.show()

        if save_fig:
            plt.savefig(os.path.join(fig_path, "significant_results.svg"), bbox_inches="tight", pad_inches=0, dpi=500)
        return

    logger.warning("Nothing was plotted. Check input data!")
    return
