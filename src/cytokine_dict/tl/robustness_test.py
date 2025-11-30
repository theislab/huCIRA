import warnings

import numpy as np
import pandas as pd
from tqdm.auto import tqdm


def _check_robustness_fractions(
    df_pivot,
    threshold_qval=0.1,  # adjusted p value
    threshold_valid=0.1,  # fraction of results required to even consider this condition. I.e. if the test only ran for one set of thresholds, then it is not very robust.
    threshold_below_alpha=0.75,  # fraction of results that need to be significant
):
    n_total = np.prod(df_pivot.shape)
    n_valid = n_total - df_pivot.isna().sum().sum()
    n_below_alpha = (
        (df_pivot < threshold_qval).sum().sum()
    )  # number of results below pval threshold, i.e., number of significant results
    frac_valid_results = n_valid / n_total
    frac_pval_below_alpha = n_below_alpha / n_valid  # fraction of significant results relative to valid results
    is_robust = (frac_pval_below_alpha > threshold_below_alpha) & (frac_valid_results > threshold_valid)
    return frac_valid_results, frac_pval_below_alpha, is_robust


def check_robustness(
    all_results,
    threshold_qval=0.1,
    threshold_valid=0.1,
    threshold_below_alpha=0.9,
):
    """Filters for robust and significant results out of original enrichments (run_enrichment_test() output)

    Returns only the enrichments that are stable across many different tests and statistically signficant.


    Parameters
    ----------
    - results
        The DataFrame output from run_enrichment_test().
    - threshold_qval
        Threshold that checks significance of results (leniently). Result is considered significant if its q-val is below this threshold.
    - threshold_valid
        The fraction of results required to even consider this condition. I.e. if the test only ran for one set of thresholds, then it is not very robust.
    - threshold_below_alpha
        The fraction of results that need to be significant


    Returns
    -------
    - robust_results
        DataFrame with robust and significant cytokine enrichments (includes min and max of NES)

    """
    all_thresholds_expression = all_results.threshold_expression.sort_values(ascending=False).unique()
    all_thresholds_lfc = sorted(all_results.threshold_lfc.unique())

    df = pd.DataFrame(index=all_thresholds_expression, columns=all_thresholds_lfc)
    df.index.rename("threshold_expression", inplace=True)
    df.columns.rename("threshold_lfc", inplace=True)

    robust_results = []

    for contrast in tqdm(all_results.contrast.unique()):
        for celltype_combo in all_results.celltype_combo.unique():
            results_ct = all_results.loc[
                (all_results.celltype_combo == celltype_combo) & (all_results.contrast == contrast)
            ]
            for cytokine in results_ct.cytokine.unique():
                results_ct_cy = results_ct.loc[results_ct.cytokine == cytokine]
                df_pivot = results_ct_cy.pivot(
                    index="threshold_expression", columns="threshold_lfc", values="FDR q-val"
                )
                with warnings.catch_warnings():
                    warnings.simplefilter(action="ignore", category=FutureWarning)
                    df_combined = pd.concat([df, df_pivot])
                df_merged = df_combined.combine_first(df_pivot)
                df_merged = df_merged.loc[~df_merged.index.duplicated()]
                df_pivot = df_merged.loc[all_thresholds_expression, all_thresholds_lfc].astype(float)
                frac_valid_results, frac_pval_below_alpha, is_robust = _check_robustness_fractions(
                    df_pivot,
                    threshold_qval=threshold_qval,
                    threshold_valid=threshold_valid,
                    threshold_below_alpha=threshold_below_alpha,
                )

                if is_robust:
                    robust_results.append(
                        (
                            celltype_combo,
                            contrast,
                            cytokine,
                            frac_valid_results,
                            frac_pval_below_alpha,
                            is_robust,
                            results_ct_cy.NES.min(),
                            results_ct_cy.NES.max(),
                            threshold_qval,
                            threshold_below_alpha,
                        )
                    )

    robust_results = pd.DataFrame(robust_results).rename(
        {
            0: "celltype_combo",
            1: "contrast",
            2: "cytokine",
            3: "frac_valid",
            4: "frac_significant",
            5: "is_robust",
            6: "NES_min",
            7: "NES_max",
            8: "qval_threshold",
            9: "threshold_frac_below_alpha",
        },
        axis=1,
    )
    return robust_results


def get_robust_significant_results(
    results,
    alphas=None,
    threshold_valid=0.1,
    threshold_below_alpha=0.9,
):
    """Filters for robust and signifcant results out of original enrichments (run_enrichment_test() output)

    Returns only the enrichments that are statistically significant (q-val), and stable across many different tests.
    Calls check_robustness for different qval thresholds to explore more stringent significance thresholds. Use for visualization of results (e.g. in a heatmap).

    Parameters
    ----------
    - results
        The DataFrame output from run_enrichment_test().
    - alphas
        List of thresholds (q-val) to check significance of results. Result is considered significant if its q-val is below this threshold.
    - threshold_valid
        The fraction of results required to even consider this condition. I.e. if the test only ran for one set of thresholds, then it is not very robust.
    - threshold_below_alpha
        The fraction of results that need to be significant


    Returns
    -------
    - results_pivot
        DataFrame with robust and significant cytokine enrichments (includes mean of NES)
    - df_annot
        DataFrame visualizing significance of a cytokine-celltype pair as stars. (0.1: *, 0.05: **, 0.01: ***)
    - results_robust
        DataFrame with robust and significant cytokine enrichments (includes min and max of NES)
    """
    # default significant values (matching significance stars)
    if alphas is None:
        alphas = [0.1, 0.05, 0.01]

    results_robust = []
    for alpha in alphas:
        results_robust.append(
            check_robustness(
                results,
                threshold_qval=alpha,
                threshold_valid=threshold_valid,
                threshold_below_alpha=threshold_below_alpha,
            )
        )

    results_robust = pd.concat(results_robust)
    results_robust = (
        results_robust.groupby(["contrast", "celltype_combo", "cytokine"])["qval_threshold"]
        .min()
        .to_frame()
        .reset_index()
    )

    results_mean = (
        results.fillna(0).groupby(["contrast", "celltype_combo", "cytokine"])["NES"].mean().to_frame().reset_index()
    )

    if len(results.contrast.unique()) > 1:
        results_pivot = results_mean.pivot(
            index=["contrast", "cytokine"],
            columns="celltype_combo",
            values="NES",
        )
    else:
        results_pivot = results_mean.pivot(
            index="cytokine",
            columns="celltype_combo",
            values="NES",
        )
        
    df_annot = results_pivot.copy()
    with warnings.catch_warnings():
        warnings.simplefilter(action="ignore", category=FutureWarning)
        df_annot[:] = ""
    for cytokine in df_annot.index:
        for celltype in df_annot.columns:
            qval = results_robust.loc[
                (results_robust.cytokine == cytokine) & (results_robust.celltype_combo == celltype)
            ].qval_threshold
            if len(qval) != 0:
                qval = qval.values[0]
                if qval == 0.1:
                    df_annot.loc[cytokine, celltype] = "*"
                if qval == 0.05:
                    df_annot.loc[cytokine, celltype] = "**"
                if qval == 0.01:
                    df_annot.loc[cytokine, celltype] = "***"

    return results_pivot, df_annot, results_robust
