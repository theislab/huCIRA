import logging
from collections.abc import Sequence
from typing import Literal

import numpy as np
import pandas as pd
from anndata import AnnData

logger = logging.getLogger(__name__)

_SCORE_METHODS = ("ulm", "mlm")


def _import_decoupler():
    """Import :mod:`decoupler` lazily so it stays an optional dependency."""
    try:
        import decoupler as dc
    except ImportError as e:  # pragma: no cover - depends on user environment
        raise ImportError(
            "Contrast-free cytokine activity scoring requires the 'decoupler' package. "
            'Install it with `pip install "hucira[activity]"` or `pip install decoupler`.'
        ) from e
    return dc


def get_cytokine_net(
    df: pd.DataFrame,
    celltype_signature: str,
    direction: Literal["upregulated", "downregulated", "both"] = "both",
    threshold_pval: float | None = 0.01,
    threshold_lfc: float | None = 1.0,
    var_names: Sequence[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Convert a gene-program DataFrame into a :mod:`decoupler` network for one cell type.

    The returned network is the *weighted* counterpart of the gene sets built by
    :func:`~hucira.tl.enrichment_test._get_genesets`: instead of discarding the
    effect size, the Human Cytokine Dictionary ``log_fc`` is kept as the signed
    edge weight, which is what linear-model methods such as ULM and MLM use.

    Filtering follows the same convention as the enrichment tests, i.e. for
    ``direction="downregulated"`` the *signed* ``threshold_lfc`` is used as an
    upper bound (pass e.g. ``-1.0``), and for ``direction="both"`` the absolute
    log-fold change is thresholded.

    Parameters
    ----------
    df : pandas.DataFrame
        Human Cytokine Dictionary (columns ``["log_fc", "adj_p_value", "cytokine"]``),
        CIP signatures (``["gene", "CIP", "celltype"]``), or a custom gene-program
        DataFrame (``["gene", "query_program", "celltype"]``).  CIP and custom
        signatures are unweighted and therefore get a weight of ``1.0``.
    celltype_signature : str
        Cell-type name matching ``df.celltype``.  Only rows of this cell type
        enter the network, which is what makes the resulting activities
        cell-type specific.
    direction : str
        One of ``"upregulated"``, ``"downregulated"`` or ``"both"``.  For
        weighted scoring ``"both"`` is usually what you want, because the
        signed weights already encode the direction.
    threshold_pval : float or None
        Maximum adjusted p-value.  Only relevant for the Human Cytokine Dictionary.
    threshold_lfc : float or None
        Log-fold change threshold.  Only relevant for the Human Cytokine Dictionary.
    var_names : sequence of str or None
        Genes present in the query data.  If given, the network is restricted to
        these genes and the overlap statistics are computed against them.

    Returns
    -------
    net : pandas.DataFrame
        Long-format network with columns ``["source", "target", "weight"]``,
        ready to be passed to :func:`decoupler.mt.ulm` or :func:`decoupler.mt.mlm`.
    net_stats : pandas.DataFrame
        Per-source gene-overlap statistics between the query data and the gene
        program, with the same column names as the enrichment tests use.
    """
    required_for_hcd = ["log_fc", "adj_p_value", "cytokine"]
    required_for_CIP = ["gene", "CIP", "celltype"]

    if set(required_for_hcd).issubset(df.columns):
        program_column = "cytokine"
        logger.info("Building weighted cytokine network of Human Cytokine Dictionary for %s.", celltype_signature)
        select = (df.adj_p_value <= threshold_pval) & (df.celltype == celltype_signature)
        if direction == "upregulated":
            select = select & (df.log_fc >= threshold_lfc)
        elif direction == "downregulated":
            select = select & (df.log_fc <= threshold_lfc)
        elif direction == "both":
            select = select & (df.log_fc.abs() >= threshold_lfc)
        else:
            raise ValueError(f"Invalid direction: {direction}.")
        net = df.loc[select, ["cytokine", "gene", "log_fc"]].copy()
        net.columns = ["source", "target", "weight"]

    elif set(required_for_CIP).issubset(df.columns):
        program_column = "CIP"
        logger.info("Building network of Cytokine-induced gene programs for %s.", celltype_signature)
        net = df.loc[df.celltype == celltype_signature, ["CIP", "gene"]].copy()
        net.columns = ["source", "target"]
        net["weight"] = 1.0

    elif "query_program" in df.columns:
        program_column = "query_program"
        logger.info("Building network of user-defined gene programs for %s.", celltype_signature)
        net = df.loc[df.celltype == celltype_signature, ["query_program", "gene"]].copy()
        net.columns = ["source", "target"]
        net["weight"] = 1.0

    else:
        raise ValueError(
            "invalid input for df parameter. You can use either the Human Cytokine Dictionary with "
            "load_human_cytokine_dict(), or our CIP signatures with load_CIP_signatures(). If you want to compute "
            "activities of custom gene sets, df must have columns: ['gene', 'query_program', 'celltype']."
        )

    net = net.dropna(subset=["source", "target", "weight"])
    # A gene may appear twice for the same source (e.g. several contrasts); keep the strongest effect.
    net = net.sort_values("weight", key=lambda weights: weights.abs(), ascending=False)
    net = net.drop_duplicates(subset=["source", "target"], keep="first")

    n_total = net.groupby("source", observed=True).size().rename("num_genes_signature")
    if var_names is not None:
        shared = net.target.isin(set(var_names))
        n_shared = net.loc[shared].groupby("source", observed=True).size().rename("num_shared_genes_signature")
        net = net.loc[shared]
    else:
        n_shared = n_total.rename("num_shared_genes_signature")

    net_stats = pd.concat([n_total, n_shared], axis=1).fillna({"num_shared_genes_signature": 0})
    net_stats["frac_shared_genes_signature"] = net_stats.num_shared_genes_signature / net_stats.num_genes_signature
    net_stats = net_stats.rename_axis(program_column).reset_index()

    net = net.sort_values(["source", "target"]).reset_index(drop=True)
    return net, net_stats


def _to_dense(X) -> np.ndarray:
    return np.asarray(X.todense() if hasattr(X, "todense") else X, dtype=np.float32)


def _aggregate(adata: AnnData, groupby: list[str], min_cells: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Average expression per group of cells (mean of the normalised matrix)."""
    obs = adata.obs
    keys = obs[groupby].astype(str)
    group_labels = keys.agg("|".join, axis=1) if len(groupby) > 1 else keys.iloc[:, 0]

    profiles, meta = {}, {}
    for group, idx in group_labels.groupby(group_labels).groups.items():
        positions = obs.index.get_indexer(idx)
        n_cells = len(positions)
        if n_cells < min_cells:
            logger.warning("Skipping group '%s': only %d cells (min_cells=%d).", group, n_cells, min_cells)
            continue
        profiles[group] = np.asarray(adata.X[positions].mean(axis=0), dtype=np.float32).ravel()
        first = keys.loc[idx[0]]
        meta[group] = {**{k: first[k] for k in groupby}, "num_cells": n_cells}

    if not profiles:
        return pd.DataFrame(columns=adata.var_names), pd.DataFrame()

    mat = pd.DataFrame.from_dict(profiles, orient="index", columns=list(adata.var_names))
    return mat, pd.DataFrame.from_dict(meta, orient="index")


def _scale(mat: pd.DataFrame) -> pd.DataFrame:
    """Z-score every gene across observations, so scores are relative within the data."""
    values = mat.to_numpy(dtype=np.float32, copy=True)
    mean = values.mean(axis=0)
    std = values.std(axis=0)
    std[std == 0] = 1.0  # constant genes become 0 and cannot drive any score
    return pd.DataFrame((values - mean) / std, index=mat.index, columns=mat.columns)


def run_one_ulm_test(
    adata: AnnData,
    df: pd.DataFrame,
    celltype_combo: tuple[str, str],
    celltype_column: str,
    groupby: str | Sequence[str] | None = None,
    direction: Literal["upregulated", "downregulated", "both"] = "both",
    # Filtering parameters for network construction
    threshold_lfc: float = 1.0,
    threshold_expression: float = 0.0,
    threshold_pval: float = 0.01,
    # Scoring parameters
    method: Literal["ulm", "mlm"] = "ulm",
    scale: bool = True,
    tmin: int = 10,
    min_cells: int = 10,
    key_added: str | None = None,
    verbose: bool = False,
    **kwargs,
) -> pd.DataFrame | None:
    """Score cytokine activity in one cell type **without a contrast**.

    Unlike :func:`run_one_enrichment_test`, which needs two conditions to rank
    genes, this estimates an activity per observation (per cell, or per group of
    cells when *groupby* is given) with a linear model, in the spirit of
    PROGENy pathway activities:

    1. "Looks up" the query cell type in the gene program (e.g. the Human
       Cytokine Dictionary) and builds a *weighted*, signed network from it
       (:func:`get_cytokine_net`).
    2. Optionally aggregates cells into groups (e.g. per donor) and z-scores
       each gene across observations.
    3. Fits, for every observation, a linear model of expression on the
       cytokine's gene weights and takes the t-value of the slope as the
       activity score (see :func:`decoupler.mt.ulm`).

    Because the network is built from the dictionary rows of *one* cell type and
    only cells annotated as the matching query cell type are scored, activities
    are cell-type specific by construction.

    Parameters
    ----------
    adata : ~anndata.AnnData
        The query AnnData object.  ``adata.X`` is expected to be normalised and
        log-transformed.
    df : pandas.DataFrame
        Human Cytokine Dictionary, CIP signatures, or custom gene-program DataFrame.
    celltype_combo : tuple of (str, str)
        Cell-type name in the query *adata* (first) and the corresponding
        cell-type name in *df* (second).
    celltype_column : str
        Column name in ``adata.obs`` that stores the cell types.
    groupby : str or sequence of str or None
        Key(s) in ``adata.obs`` to aggregate cells by before scoring, e.g.
        ``"donor_id"`` or ``["donor_id", "timepoint"]``.  Strongly recommended:
        single-cell profiles are sparse, so per-cell scores are noisy.  If
        ``None``, every cell is scored individually.
    direction : str
        One of ``"upregulated"``, ``"downregulated"`` or ``"both"``.  Defaults to
        ``"both"``, since the signed weights already carry the direction.
    threshold_lfc : float
        Log-fold change threshold used to build the network.
    threshold_expression : float
        Minimum mean gene expression across the cells of this cell type.  Note
        that this changes the gene universe of the regression and therefore the
        scores.
    threshold_pval : float
        Maximum adjusted p-value used to build the network.
    method : str
        ``"ulm"`` (univariate, one cytokine at a time) or ``"mlm"`` (multivariate,
        all cytokines jointly).  ``"mlm"`` accounts for the strong target overlap
        between related cytokines but becomes unstable when signatures are nearly
        collinear -- inspect :func:`decoupler.pp.net_corr` before trusting it.
    scale : bool
        Whether to z-score each gene across observations before scoring.  With
        ``scale=False`` the scores are dominated by baseline gene abundance; with
        ``scale=True`` they express deviation relative to the rest of the data
        and are therefore only comparable *within* this call.
    tmin : int
        Minimum number of shared target genes for a cytokine to be scored.
        Cytokines below this are dropped (they will be missing from the output).
    min_cells : int
        Minimum number of cells per group when *groupby* is used.
    key_added : str or None
        If given and *groupby* is ``None``, scores and adjusted p-values are also
        written to ``adata.obsm[f"{key_added}_score"]`` and
        ``adata.obsm[f"{key_added}_padj"]`` of the *input* object, with ``NaN``
        for cells of other cell types.
    verbose : bool
        Whether :mod:`decoupler` should print progress messages.
    kwargs
        Additional keyword arguments forwarded to the :mod:`decoupler` method.

    Returns
    -------
    pandas.DataFrame
        Long-format DataFrame with one row per observation and cytokine, holding
        ``score`` (the t-value) and ``padj``, plus the parameters used.  Returns
        ``None`` if the cell type is absent or nothing could be scored.
    """
    if method not in _SCORE_METHODS:
        raise ValueError(f"Invalid method: {method}. Choose one of {_SCORE_METHODS}.")

    dc = _import_decoupler()

    celltype_adata, celltype_signature = celltype_combo

    # allows potential loop of celltype combos to continue
    if celltype_adata not in adata.obs[celltype_column].unique():
        logger.warning(
            "'%s' is not present in celltype_column (%s) of query adata. Skipping activity test of this celltype.",
            celltype_adata,
            celltype_column,
        )
        return None

    adata_input = adata
    # filter for cell type -- cells are only ever scored with their own cell type's network
    adata = adata[adata.obs[celltype_column] == celltype_adata]
    # filter based on gene expression
    adata = adata[:, np.asarray(adata.X.mean(axis=0)).ravel() >= threshold_expression]

    net, net_stats = get_cytokine_net(
        df=df,
        celltype_signature=celltype_signature,
        direction=direction,
        threshold_pval=threshold_pval,
        threshold_lfc=threshold_lfc,
        var_names=adata.var_names,
    )
    program_column = net_stats.columns[0]
    if net.empty:
        logger.warning(
            "No network edges left for celltype '%s' (threshold_lfc=%s, threshold_pval=%s). Skipping.",
            celltype_signature,
            threshold_lfc,
            threshold_pval,
        )
        return None

    # decoupler raises if no source reaches tmin, so drop small signatures here and
    # skip the cell type gracefully -- this keeps threshold sweeps running.
    n_shared = net.groupby("source", observed=True).size()
    net = net.loc[net.source.isin(n_shared[n_shared >= tmin].index)]
    if net.empty:
        logger.warning(
            "No %s of '%s' reached tmin=%d shared target genes. Skipping.", program_column, celltype_signature, tmin
        )
        return None

    if groupby is None:
        mat = pd.DataFrame(_to_dense(adata.X), index=adata.obs_names, columns=adata.var_names)
        obs_meta = pd.DataFrame({"num_cells": 1}, index=mat.index)
        groupby_cols: list[str] = []
    else:
        groupby_cols = [groupby] if isinstance(groupby, str) else list(groupby)
        missing = [key for key in groupby_cols if key not in adata.obs.columns]
        if missing:
            raise KeyError(f"groupby key(s) not found in adata.obs: {missing}")
        mat, obs_meta = _aggregate(adata, groupby_cols, min_cells=min_cells)
        if mat.empty:
            logger.warning("No group of '%s' passed min_cells=%d. Skipping.", celltype_adata, min_cells)
            return None

    if mat.shape[0] < 2 and scale:
        logger.warning(
            "Only %d observation(s) for '%s'; z-scoring genes is not meaningful. Continuing with scale=False.",
            mat.shape[0],
            celltype_adata,
        )
        scale = False
    if scale:
        mat = _scale(mat)

    logger.info(
        "Running %s for '%s' (%s): %d observations, %d genes, %d cytokines.",
        method,
        celltype_adata,
        celltype_signature,
        mat.shape[0],
        mat.shape[1],
        net.source.nunique(),
    )
    score, padj = getattr(dc.mt, method)(mat, net, tmin=tmin, verbose=verbose, **kwargs)

    results = (
        score.rename_axis("observation")
        .reset_index()
        .melt(id_vars="observation", var_name=program_column, value_name="score")
    )
    results["padj"] = (
        padj.rename_axis("observation")
        .reset_index()
        .melt(id_vars="observation", var_name=program_column, value_name="padj")["padj"]
    )

    results = results.merge(obs_meta.rename_axis("observation").reset_index(), on="observation", how="left")
    results["celltype_adata"] = celltype_adata
    results["celltype_signature"] = celltype_signature
    results["celltype_combo"] = f"{celltype_adata} ({celltype_signature})"
    results["method"] = method
    results["direction"] = direction
    results["scale"] = scale
    results["threshold_pval"] = threshold_pval
    results["threshold_lfc"] = threshold_lfc
    results["threshold_expression"] = threshold_expression
    results["tmin"] = tmin
    results["num_genes_universe"] = mat.shape[1]
    # deliberately not called "groupby": that would shadow DataFrame.groupby on the result
    results["groupby_keys"] = "|".join(groupby_cols) if groupby_cols else None
    results = results.merge(net_stats, on=program_column, how="left")

    if key_added is not None:
        if groupby is not None:
            logger.warning("key_added is ignored when groupby is used, because scores are not per cell.")
        else:
            for name, values in (("score", score), ("padj", padj)):
                obsm_key = f"{key_added}_{name}"
                stored = adata_input.obsm.get(obsm_key)
                if not isinstance(stored, pd.DataFrame):
                    stored = pd.DataFrame(index=adata_input.obs_names, dtype=np.float32)
                for column in values.columns:
                    if column not in stored.columns:
                        stored[column] = np.nan
                stored.loc[values.index, values.columns] = values.to_numpy()
                adata_input.obsm[obsm_key] = stored

    return results


def run_all_ulm_test(
    adata: AnnData,
    df: pd.DataFrame,
    celltype_combos: list[tuple[str, str]],
    celltype_column: str,
    groupby: str | Sequence[str] | None = None,
    direction: Literal["upregulated", "downregulated", "both"] = "both",
    # Filtering parameters for network construction
    threshold_lfc: float | Sequence[float] = 1.0,
    threshold_expression: float | Sequence[float] = 0.0,
    threshold_pval: float = 0.01,
    # Scoring parameters
    method: Literal["ulm", "mlm"] = "ulm",
    scale: bool = True,
    tmin: int = 10,
    min_cells: int = 10,
    key_added: str | None = None,
    verbose: bool = False,
    **kwargs,
) -> pd.DataFrame:
    """Score contrast-free cytokine activity across cell types and thresholds.

    Wrapper around :func:`run_one_ulm_test` that loops over cell types and over
    several threshold values, so the sensitivity of an activity to the network
    definition can be assessed the same way :func:`run_all_enrichment_test`
    does for enrichment.

    Every cell type is scored against its own slice of the gene program, and the
    results are concatenated into one long-format DataFrame.

    Parameters
    ----------
    adata : ~anndata.AnnData
        The query AnnData object.
    df : pandas.DataFrame
        Human Cytokine Dictionary, CIP signatures, or custom gene-program DataFrame.
    celltype_combos : list of tuple
        List of ``(query_celltype, df_celltype)`` pairs, e.g. from
        :func:`create_celltype_combos`.
    celltype_column : str
        Column name in ``adata.obs`` that stores the cell types.
    groupby : str or sequence of str or None
        Key(s) in ``adata.obs`` to aggregate cells by before scoring.
    direction : str
        One of ``"upregulated"``, ``"downregulated"`` or ``"both"``.
    threshold_lfc : float or sequence of float
        Log-fold change threshold(s) for building the network.
    threshold_expression : float or sequence of float
        Mean gene-expression threshold(s).
    threshold_pval : float
        Maximum adjusted p-value for the network.
    method : str
        ``"ulm"`` or ``"mlm"``.
    scale : bool
        Whether to z-score each gene across observations before scoring.
    tmin : int
        Minimum number of shared target genes for a cytokine to be scored.
    min_cells : int
        Minimum number of cells per group when *groupby* is used.
    key_added : str or None
        Forwarded to :func:`run_one_ulm_test`; only sensible for a single
        threshold combination.
    verbose : bool
        Whether :mod:`decoupler` should print progress messages.
    kwargs
        Additional keyword arguments forwarded to the :mod:`decoupler` method.

    Returns
    -------
    pandas.DataFrame
        Concatenated long-format activity scores for all cell types and thresholds.
    """
    thresholds_lfc = [threshold_lfc] if np.isscalar(threshold_lfc) else list(threshold_lfc)
    thresholds_expression = [threshold_expression] if np.isscalar(threshold_expression) else list(threshold_expression)

    all_results = []
    for celltype_combo in celltype_combos:
        for lfc in thresholds_lfc:
            for expr in thresholds_expression:
                results = run_one_ulm_test(
                    adata=adata,
                    df=df,
                    celltype_combo=celltype_combo,
                    celltype_column=celltype_column,
                    groupby=groupby,
                    direction=direction,
                    threshold_lfc=lfc,
                    threshold_expression=expr,
                    threshold_pval=threshold_pval,
                    method=method,
                    scale=scale,
                    tmin=tmin,
                    min_cells=min_cells,
                    key_added=key_added,
                    verbose=verbose,
                    **kwargs,
                )
                if results is not None:
                    all_results.append(results)

    if not all_results:
        logger.warning("No cell type could be scored; returning an empty DataFrame.")
        return pd.DataFrame()

    return pd.concat(all_results, axis=0, ignore_index=True)
