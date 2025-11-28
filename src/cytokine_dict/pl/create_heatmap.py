import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
import seaborn as sns


def format_celltype_combo_labels(x):
    if isinstance(x, (list, np.ndarray, pd.Index)):
        return [format_celltype_combo_labels(xx) for xx in x]
    return x.split(" (")[0]


def format_program_labels(x):
    if isinstance(x, (list, np.ndarray, pd.Index)):
        return [format_program_labels(xx) for xx in x]
    return x.split(" ")[1]


def create_heatmap(
    results,
    robust_results,
    celltypes,
    cytokines,
    path = None,
    figsize = None,
):
    df_nes = results.loc[(results.cytokine.isin(cytokines) & results.celltype_combo.isin(celltypes))]\
        .groupby(["celltype_combo", "cytokine"])["NES"].mean()\
        .to_frame()\
        .reset_index()\
        .pivot(index="cytokine", columns="celltype_combo", values="NES")
    print(df_nes)
    df_nes = df_nes.loc[:, celltypes[np.isin(celltypes, df_nes.columns)]]
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=FutureWarning)
        df_annot = df_nes.copy()
        df_annot.loc[:, :] = ""
        
    for program in df_annot.index:
        for celltype in df_annot.columns:
            select = ((robust_results.cytokine == program) & (robust_results.celltype_combo == celltype))
            alpha = robust_results.loc[select, "qval_threshold"]
            if len(alpha) == 0:
                df_annot.loc[program, celltype] = ""
            else:
                assert len(alpha) == 1
                alpha = alpha.values[0]
                if alpha == 0.1:
                    df_annot.loc[program, celltype] = "*"
                elif alpha == 0.05:
                    df_annot.loc[program, celltype] = "**"
                elif alpha == 0.01:
                    df_annot.loc[program, celltype] = "***"
    
    if figsize is None:
        fig, ax = plt.subplots(1, 1, figsize=(0.75, 0.75))
    else:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.set_facecolor("lightgray")

    print(df_nes)
    print(df_annot)
    print(print(df_nes.dtypes))
    df_nes = df_nes.apply(pd.to_numeric, errors='coerce')
    heatmap_out = sns.heatmap(
        df_nes,
        center=0,
        #vmin=-df_blood.NES.max(),
        #vmax=df_blood.NES.max(),
        cmap="RdBu_r",
        square=True,
        annot=df_annot,
        annot_kws={"fontsize": 6, "family": "sans-serif"},
        fmt="",
        ax=ax,
        cbar=False,
        cbar_kws={"shrink": 0.5},
        lw=0.5,
        linecolor="white",
    )
    plt.xlabel("")
    plt.ylabel("")
    ax.tick_params(axis='both', labelright=False, labelleft=True, length=0)
    ax.set_xticks(0.5+np.arange(df_nes.shape[1]), format_celltype_combo_labels(df_nes.columns), family="sans-serif", fontsize=6, rotation=90)
    ax.set_yticks(0.5+np.arange(df_nes.shape[0]), format_program_labels(df_nes.index), family="sans-serif", fontsize=6)
   
    if path:
        plt.savefig(path, dpi=400, bbox_inches="tight", pad_inches=0)
    plt.show()
    return df_nes, df_annot
