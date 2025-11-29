from tqdm.auto import tqdm
import warnings
import numpy as np
import pandas as pd
import os



def format_cytokine_names(x):
    if isinstance(x, (list, np.ndarray, pd.Index)):
        return [format_cytokine_names(_x) for _x in x]
    text = x.get_text() if hasattr(x, 'get_text') else x
    text = text.replace("beta", r"$\beta$")
    text = text.replace("alpha", r"$\alpha$")
    text = text.replace("gamma", r"$\gamma$")
    text = text.replace("lambda", r"$\lambda$")
    text = text.replace("omega", r"$\omega$")
    return text


'''
def format_celltype_names(x):
    if isinstance(x, (list, np.ndarray, pd.Index)):
        return [format_celltype_names(_x) for _x in x]
    text = x.get_text() if hasattr(x, 'get_text') else x
    # text = text.replace("CD4-positive, alpha-beta T cell", r"CD4, $\alpha$-$\beta$ T")
    text = text.replace("CD4-positive, alpha-beta T cell", r"CD4 T")
    # text = text.replace("CD8-positive, alpha-beta T cell", r"CD8, $\alpha$-$\beta$ T")
    text = text.replace("CD8-positive, alpha-beta T cell", r"CD8 T")
    text = text.replace("conventional dendritic cell", "cDC")
    text = text.replace("natural killer cell", "NK")
    if not text.startswith("NK"):
        text = text.split("(")[0]
    else:
        #text = text.replace("CD56-bright_NK", "CD56hi")
        #text = text.replace("CD56-dim_NK", "CD56low")
        #text = text.replace("NK_CD56hi", "CD56hi")
        text = text.replace("NK_CD56hi", "hi")
        #text = text.replace("NK_CD56low", "CD56low")
        text = text.replace("NK_CD56low", "low")
    if text.startswith("non-classical monocyte"):
        text =  "nc. Mono."
    elif text.startswith("classical monocyte"):
        text =  "c. Mono."
    return text
'''


def plot_significant_results(results_pivot, df_annot, celltype_ordered, fontsize=6, save_fig=False, fig_path=""):

    # Why pivot at all?
    results_pivot = results_pivot.T.loc[CELLTYPE_ORDER].T

    # Plot
    fig, ax = plt.subplots(1,1,figsize=(8, 8)) #?
    sns.heatmap(
        results_pivot,
        square=True,
        annot=df_annot,
        cmap="RdBu_r",
        center=0,
        # vmin=vmin,
        # vmax=vmax,
        annot_kws={"fontsize": 6, "family": "sans-serif"},
        linecolor="white",
        fmt="",
        linewidths=0.5,
        cbar=False,
        ax=ax,
    )

    plt.xlabel("")
    plt.ylabel("")
    ax.set_facecolor("lightgray")
    ax.tick_params(axis='both', which='both', length=0)

    # ToDo: Decide either delete function or make more general.
    # celltype_labels = format_celltype_names(results_pivot.columns)
    celltype_labels = results_pivot.columns
    plt.xticks(0.5+np.arange(results_pivot.shape[1]), celltype_labels, fontsize=fontsize, family="sans-serif", rotation=90, ha="center")
    
    cytokine_labels = format_cytokine_names(results_pivot.index)
    plt.yticks(0.5+np.arange(results_pivot.shape[0]), cytokine_labels, fontsize=fontsize, family="sans-serif", rotation=0, ha="right")
    plt.show()

    if save_fig:
        plt.savefig(
            os.path.join(fig_path, "significant_results.svg"), 
            bbox_inches="tight", 
            pad_inches=0, 
            dpi=500,
            # transparent=True,
        )

    return

