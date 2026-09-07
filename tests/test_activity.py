import anndata as ad
import numpy as np
import pandas as pd
import pytest

import hucira as hc

N_GENES = 60
CELLTYPES = {"B cell": "B", "T cell": "T"}
# each cell type gets its own, non-overlapping signature for the same cytokine
SIGNATURES = {"B": list(range(0, 15)), "T": list(range(30, 45))}


@pytest.fixture
def hcd():
    """Minimal Human Cytokine Dictionary: IL-TEST is cell-type specific, IL-OTHER is noise."""
    rows = []
    for celltype, genes in SIGNATURES.items():
        for gene_i in genes:
            rows.append(("IL-TEST", celltype, f"gene_{gene_i}", 2.0, 0.001))
        for gene_i in range(45, 60):
            rows.append(("IL-OTHER", celltype, f"gene_{gene_i}", 1.5, 0.001))
        # below the p-value / lfc thresholds, must be filtered out
        rows.append(("IL-TEST", celltype, "gene_59", 3.0, 0.5))
        rows.append(("IL-TEST", celltype, "gene_58", 0.1, 0.001))
    return pd.DataFrame(rows, columns=["cytokine", "celltype", "gene", "log_fc", "adj_p_value"])


@pytest.fixture
def adata_activity():
    """Two cell types x two donors; only donor_1 has the cell-type-specific program on."""
    rng = np.random.default_rng(0)
    obs = pd.DataFrame(
        [(celltype, donor, i) for celltype in CELLTYPES for donor in ("donor_0", "donor_1") for i in range(30)],
        columns=["celltype", "donor_id", "rep"],
    )
    obs.index = [f"cell_{i}" for i in range(len(obs))]
    X = rng.normal(loc=5.0, scale=1.0, size=(len(obs), N_GENES)).astype(np.float32)

    on = (obs.donor_id == "donor_1").to_numpy()
    for celltype, signature_name in CELLTYPES.items():
        rows = on & (obs.celltype == celltype).to_numpy()
        X[np.ix_(rows, SIGNATURES[signature_name])] += 4.0

    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"gene_{i}" for i in range(N_GENES)]))


def test_get_cytokine_net_filters_and_keeps_signed_weights(hcd, adata_activity):
    net, net_stats = hc.get_cytokine_net(
        hcd, celltype_signature="B", threshold_lfc=1.0, threshold_pval=0.01, var_names=adata_activity.var_names
    )

    il_test = net.loc[net.source == "IL-TEST"]
    assert set(il_test.target) == {f"gene_{i}" for i in SIGNATURES["B"]}  # T-cell genes excluded
    assert (il_test.weight == 2.0).all()
    assert list(net.columns) == ["source", "target", "weight"]
    assert net_stats.set_index("cytokine").loc["IL-TEST", "frac_shared_genes_signature"] == 1.0


def test_get_cytokine_net_direction_and_missing_genes(hcd):
    hcd = pd.concat([hcd, pd.DataFrame([("IL-TEST", "B", "gene_absent", -2.0, 0.001)], columns=hcd.columns)])

    up, _ = hc.get_cytokine_net(hcd, "B", direction="upregulated", threshold_lfc=1.0)
    down, _ = hc.get_cytokine_net(hcd, "B", direction="downregulated", threshold_lfc=-1.0)
    assert (up.weight > 0).all()
    assert list(down.target) == ["gene_absent"]

    both, stats = hc.get_cytokine_net(hcd, "B", direction="both", threshold_lfc=1.0, var_names=["gene_0", "gene_1"])
    assert set(both.target) == {"gene_0", "gene_1"}  # restricted to var_names
    assert stats.set_index("cytokine").loc["IL-TEST", "num_shared_genes_signature"] == 2

    with pytest.raises(ValueError, match="Invalid direction"):
        hc.get_cytokine_net(hcd, "B", direction="sideways")


def test_run_one_ulm_test_is_celltype_specific(hcd, adata_activity):
    results = hc.run_all_ulm_test(
        adata=adata_activity,
        df=hcd,
        celltype_combos=list(CELLTYPES.items()),
        celltype_column="celltype",
        groupby="donor_id",
        tmin=5,
        min_cells=5,
    )

    assert set(results.celltype_combo) == {"B cell (B)", "T cell (T)"}
    assert set(results.cytokine) == {"IL-TEST", "IL-OTHER"}
    assert (results.num_cells == 30).all()
    assert results.groupby(["celltype_combo", "cytokine", "observation"]).size().max() == 1

    scores = results.pivot_table(index=["celltype_combo", "cytokine"], columns="observation", values="score")

    # donor_1 has the program on in both cell types, and it is picked up in each
    for celltype_combo in ("B cell (B)", "T cell (T)"):
        assert scores.loc[(celltype_combo, "IL-TEST"), "donor_1"] > 0
        assert scores.loc[(celltype_combo, "IL-TEST"), "donor_0"] < 0

    # the B-cell IL-TEST signature must not be used to score T cells: swapping the
    # dictionary cell type destroys the signal, which is what "cell-type specific" means
    swapped = hc.run_one_ulm_test(
        adata=adata_activity,
        df=hcd,
        celltype_combo=("T cell", "B"),  # T cells scored with the B-cell signature
        celltype_column="celltype",
        groupby="donor_id",
        tmin=5,
        min_cells=5,
    )
    swapped_score = swapped.set_index(["cytokine", "observation"]).loc[("IL-TEST", "donor_1"), "score"]
    assert swapped_score < scores.loc[("T cell (T)", "IL-TEST"), "donor_1"]


def test_run_one_ulm_test_writes_obsm_with_nan_outside_celltype(hcd, adata_activity):
    hc.run_one_ulm_test(
        adata=adata_activity,
        df=hcd,
        celltype_combo=("B cell", "B"),
        celltype_column="celltype",
        tmin=5,
        key_added="cytokine_activity",
    )

    scores = adata_activity.obsm["cytokine_activity_score"]
    is_b_cell = (adata_activity.obs.celltype == "B cell").to_numpy()
    assert scores.shape == (adata_activity.n_obs, 2)
    assert scores.loc[is_b_cell].notna().all().all()
    assert scores.loc[~is_b_cell].isna().all().all()
    assert adata_activity.obsm["cytokine_activity_padj"].loc[is_b_cell].notna().all().all()

    # a second cell type fills in the remaining rows of the same obsm frame
    hc.run_one_ulm_test(
        adata=adata_activity,
        df=hcd,
        celltype_combo=("T cell", "T"),
        celltype_column="celltype",
        tmin=5,
        key_added="cytokine_activity",
    )
    assert adata_activity.obsm["cytokine_activity_score"].notna().all().all()


def test_run_one_ulm_test_multiple_groupby_keys_and_mlm(hcd, adata_activity):
    results = hc.run_one_ulm_test(
        adata=adata_activity,
        df=hcd,
        celltype_combo=("B cell", "B"),
        celltype_column="celltype",
        groupby=["donor_id", "celltype"],
        method="mlm",
        tmin=5,
        min_cells=5,
    )

    assert set(results.observation) == {"donor_0|B cell", "donor_1|B cell"}
    assert set(results.donor_id) == {"donor_0", "donor_1"}
    assert (results.method == "mlm").all()
    assert (results.groupby_keys == "donor_id|celltype").all()


def test_run_one_ulm_test_skips_and_validates(hcd, adata_activity):
    assert (
        hc.run_one_ulm_test(
            adata=adata_activity,
            df=hcd,
            celltype_combo=("Missing cell", "B"),
            celltype_column="celltype",
        )
        is None
    )
    # no cytokine reaches tmin
    assert (
        hc.run_one_ulm_test(
            adata=adata_activity, df=hcd, celltype_combo=("B cell", "B"), celltype_column="celltype", tmin=500
        )
        is None
    )
    # no network edge survives the threshold
    assert (
        hc.run_one_ulm_test(
            adata=adata_activity,
            df=hcd,
            celltype_combo=("B cell", "B"),
            celltype_column="celltype",
            threshold_lfc=10.0,
        )
        is None
    )
    with pytest.raises(ValueError, match="Invalid method"):
        hc.run_one_ulm_test(
            adata=adata_activity, df=hcd, celltype_combo=("B cell", "B"), celltype_column="celltype", method="gsea"
        )
    with pytest.raises(KeyError, match="groupby"):
        hc.run_one_ulm_test(
            adata=adata_activity,
            df=hcd,
            celltype_combo=("B cell", "B"),
            celltype_column="celltype",
            groupby="not_a_column",
        )


def test_threshold_sweep_and_sparse_input(hcd, adata_activity):
    from scipy.sparse import csr_matrix

    adata_activity.X = csr_matrix(adata_activity.X)
    results = hc.run_all_ulm_test(
        adata=adata_activity,
        df=hcd,
        celltype_combos=list(CELLTYPES.items()),
        celltype_column="celltype",
        groupby="donor_id",
        threshold_lfc=[1.0, 1.8],
        threshold_expression=[0.0, 1.0],
        tmin=5,
        min_cells=5,
    )

    assert set(results.threshold_lfc) == {1.0, 1.8}
    assert set(results.threshold_expression) == {0.0, 1.0}
    # IL-OTHER (lfc 1.5) drops out of the stricter network, IL-TEST (lfc 2.0) stays
    strict = results.loc[results.threshold_lfc == 1.8]
    assert set(strict.cytokine) == {"IL-TEST"}
