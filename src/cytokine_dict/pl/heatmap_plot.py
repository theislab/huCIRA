from tqdm.auto import tqdm
import warnings
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

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



def plot_significant_results(robust_results_dict=None, 
                             results_pivot=None, 
                             df_annot=None, 
                             selected_celltypes=None, 
                             selected_cytokines=None, 
                             fontsize=6, 
                             save_fig=False, 
                             fig_path=""):

    """ Optional heatmap plotting aid: Plots either the robust results from a dict of contrasts or individually per contrast.
    
    Parameters
    ----------
    - robust_results_dict:
        robust enrichment score dictionary from get_significant_results(). If this argument is present it has precedence over results_pivot and df_annot.
    - results_pivot:
        pandas DataFrame of robust enrichment for results from one contrast  
    - df_annot:
        pandas DataFrame of robust enrichment significance annotations for results from one contrast
    - selected_celltypes:
        Can choose to only visualize selected celltypes out of available from robust results. Must be in robust results, otherwise error.
    - selected_cytokines:
        Can choose to only visualize selected celltypes out of available from robust results. Must be in robust results, otherwise error.
    
    Returns
    ----------
    - Nothing. Plotting function only

    """
    
    
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
                ax=ax
            )

            ax.set_title(contrast, fontsize=10)
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.set_facecolor("lightgray")
            ax.tick_params(axis='both', which='both', length=0)

            # Axis labels
            ax.set_xticks(0.5 + np.arange(pivot.shape[1]))
            ax.set_xticklabels(pivot.columns, fontsize=fontsize, rotation=90, ha="center")
            ax.set_yticks(0.5 + np.arange(pivot.shape[0]))
            ax.set_yticklabels(format_cytokine_names(pivot.index), fontsize=fontsize, rotation=0, ha="right")
                    
        if save_fig:
            plt.savefig(
                os.path.join(fig_path, f"all_contrasts_significant_results.svg"),
                bbox_inches="tight",
                pad_inches=0,
                dpi=500
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
        
        fig, ax = plt.subplots(1, 1)
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
            ax=ax
        )
        ax.set_title("Contrast1_vs_Contrast2", fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_facecolor("lightgray")
        ax.tick_params(axis='both', which='both', length=0)

        # Axis labels
        ax.set_xticks(0.5 + np.arange(results_pivot.shape[1]))
        ax.set_xticklabels(results_pivot.columns, fontsize=fontsize, rotation=90, ha="center")
        ax.set_yticks(0.5 + np.arange(results_pivot.shape[0]))
        ax.set_yticklabels(format_cytokine_names(results_pivot.index), fontsize=fontsize, rotation=0, ha="right")

        plt.show()

        if save_fig:
            plt.savefig(
                os.path.join(fig_path, "significant_results.svg"),
                bbox_inches="tight",
                pad_inches=0,
                dpi=500
            )
        return

    print("Nothing was plotted. Check input data!")
    return

