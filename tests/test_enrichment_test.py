"""Regression tests for hucira.tl.enrichment_test."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from hucira.tl.enrichment_test import (
    _compute_mu_and_sigma,
    _compute_ranking_statistic,
    _compute_s2n,
    _get_genesets,
    run_all_enrichment_test,
    run_one_enrichment_test,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

N_GENES = 100
RNG = np.random.default_rng(42)


@pytest.fixture
def gene_names():
    return [f"GENE{i}" for i in range(N_GENES)]


@pytest.fixture
def adata(gene_names):
    """Synthetic AnnData with two cell types and two disease conditions.

    Cell-type A / condition_1 has higher expression for the first 30 genes,
    so those genes should have a positive S2N for condition_1 vs condition_2
    within cell-type A.
    """
    n_cells = 200
    n_genes = len(gene_names)

    X = RNG.random((n_cells, n_genes)).astype(np.float32)
    # Boost first 30 genes in the first 50 cells (cell_type A, condition_1)
    X[:50, :30] += 3.0

    obs = pd.DataFrame(
        {
            "cell_type": (["A"] * 100) + (["B"] * 100),
            "disease_state": (["condition_1"] * 50 + ["condition_2"] * 50) * 2,
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    obs["cell_type"] = obs["cell_type"].astype("category")
    obs["disease_state"] = obs["disease_state"].astype("category")

    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=pd.DataFrame(index=gene_names),
    )
    return adata


@pytest.fixture
def hcd_df(gene_names):
    """Synthetic Human Cytokine Dictionary with two cytokines."""
    rows = []
    # Cytokine CYT1: first 20 genes upregulated in cell type A
    for g in gene_names[:20]:
        rows.append({"gene": g, "cytokine": "CYT1", "celltype": "A_ref", "log_fc": 2.0, "adj_p_value": 0.001})
    # Cytokine CYT2: genes 10-30 upregulated in cell type A
    for g in gene_names[10:30]:
        rows.append({"gene": g, "cytokine": "CYT2", "celltype": "A_ref", "log_fc": 1.5, "adj_p_value": 0.005})
    # Small cytokine that won't pass min_size filter (only 5 genes)
    for g in gene_names[80:85]:
        rows.append({"gene": g, "cytokine": "CYT_SMALL", "celltype": "A_ref", "log_fc": 2.0, "adj_p_value": 0.001})
    return pd.DataFrame(rows)


@pytest.fixture
def cip_df(gene_names):
    """Synthetic CIP signatures DataFrame."""
    rows = []
    for g in gene_names[:15]:
        rows.append({"gene": g, "CIP": "CIP1", "celltype": "A_ref"})
    for g in gene_names[15:30]:
        rows.append({"gene": g, "CIP": "CIP2", "celltype": "A_ref"})
    return pd.DataFrame(rows)


@pytest.fixture
def custom_df(gene_names):
    """Synthetic custom gene-program DataFrame."""
    rows = []
    for g in gene_names[:20]:
        rows.append({"gene": g, "query_program": "PROG1", "celltype": "A_ref"})
    for g in gene_names[20:40]:
        rows.append({"gene": g, "query_program": "PROG2", "celltype": "A_ref"})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Tests: _get_genesets
# ---------------------------------------------------------------------------


class TestGetGenesets:
    def test_hcd_returns_correct_keys(self, adata, hcd_df):
        gene_set_dict, gene_set_df = _get_genesets(
            adata, hcd_df, celltype_signature="A_ref", direction="upregulated", threshold_pval=0.01, threshold_lfc=1.0
        )
        assert set(gene_set_dict.keys()) == {"CYT1", "CYT2", "CYT_SMALL"}
        assert "cytokine" in gene_set_df.columns

    def test_hcd_gene_overlap(self, adata, hcd_df):
        gene_set_dict, gene_set_df = _get_genesets(
            adata, hcd_df, celltype_signature="A_ref", direction="upregulated", threshold_pval=0.01, threshold_lfc=1.0
        )
        # All genes in HCD exist in adata, so overlap should be complete
        assert len(gene_set_dict["CYT1"]) == 20
        assert len(gene_set_dict["CYT2"]) == 20

    def test_hcd_direction_filter(self, adata, hcd_df, gene_names):
        # Add some downregulated genes
        down_rows = [
            {"gene": g, "cytokine": "CYT_DOWN", "celltype": "A_ref", "log_fc": -2.0, "adj_p_value": 0.001}
            for g in gene_names[50:65]
        ]
        hcd_with_down = pd.concat([hcd_df, pd.DataFrame(down_rows)], ignore_index=True)

        # upregulated should NOT include CYT_DOWN
        gsd_up, _ = _get_genesets(
            adata,
            hcd_with_down,
            celltype_signature="A_ref",
            direction="upregulated",
            threshold_pval=0.01,
            threshold_lfc=1.0,
        )
        assert "CYT_DOWN" not in gsd_up

        # downregulated should include CYT_DOWN but not CYT1
        gsd_down, _ = _get_genesets(
            adata,
            hcd_with_down,
            celltype_signature="A_ref",
            direction="downregulated",
            threshold_pval=0.01,
            threshold_lfc=1.0,
        )
        assert "CYT_DOWN" in gsd_down
        assert "CYT1" not in gsd_down

        # both should include all
        gsd_both, _ = _get_genesets(
            adata, hcd_with_down, celltype_signature="A_ref", direction="both", threshold_pval=0.01, threshold_lfc=1.0
        )
        assert "CYT_DOWN" in gsd_both
        assert "CYT1" in gsd_both

    def test_hcd_invalid_direction_raises(self, adata, hcd_df):
        with pytest.raises(ValueError, match="Invalid direction"):
            _get_genesets(
                adata, hcd_df, celltype_signature="A_ref", direction="invalid", threshold_pval=0.01, threshold_lfc=1.0
            )

    def test_hcd_lfc_threshold_filters(self, adata, hcd_df):
        # With high LFC threshold, CYT2 (lfc=1.5) should be excluded
        gsd, _ = _get_genesets(
            adata, hcd_df, celltype_signature="A_ref", direction="upregulated", threshold_pval=0.01, threshold_lfc=1.8
        )
        assert "CYT1" in gsd
        assert "CYT2" not in gsd

    def test_hcd_pval_threshold_filters(self, adata, hcd_df):
        # With very strict pval threshold, only CYT1 (pval=0.001) passes
        gsd, _ = _get_genesets(
            adata, hcd_df, celltype_signature="A_ref", direction="upregulated", threshold_pval=0.002, threshold_lfc=1.0
        )
        assert "CYT1" in gsd
        assert "CYT2" not in gsd

    def test_cip_returns_correct_keys(self, adata, cip_df):
        gene_set_dict, gene_set_df = _get_genesets(adata, cip_df, celltype_signature="A_ref")
        assert set(gene_set_dict.keys()) == {"CIP1", "CIP2"}
        assert "CIP" in gene_set_df.columns

    def test_cip_gene_counts(self, adata, cip_df):
        gene_set_dict, _ = _get_genesets(adata, cip_df, celltype_signature="A_ref")
        assert len(gene_set_dict["CIP1"]) == 15
        assert len(gene_set_dict["CIP2"]) == 15

    def test_custom_returns_correct_keys(self, adata, custom_df):
        gene_set_dict, gene_set_df = _get_genesets(adata, custom_df, celltype_signature="A_ref")
        assert set(gene_set_dict.keys()) == {"PROG1", "PROG2"}
        assert "query_program" in gene_set_df.columns

    def test_invalid_df_raises(self, adata):
        bad_df = pd.DataFrame({"col1": [1], "col2": [2]})
        with pytest.raises(ValueError, match="invalid input"):
            _get_genesets(adata, bad_df, celltype_signature="A_ref")

    def test_gene_set_df_has_overlap_info(self, adata, hcd_df):
        _, gene_set_df = _get_genesets(
            adata, hcd_df, celltype_signature="A_ref", direction="upregulated", threshold_pval=0.01, threshold_lfc=1.0
        )
        expected_cols = {"cytokine", "num_genes_signature", "num_shared_genes_signature", "frac_shared_genes_signature"}
        assert expected_cols.issubset(gene_set_df.columns)
        # All genes exist in adata so frac should be 1.0
        assert (gene_set_df["frac_shared_genes_signature"] == 1.0).all()

    def test_partial_gene_overlap(self, adata, gene_names):
        """When the HCD contains genes absent from adata, overlap should be partial."""
        rows = []
        for g in gene_names[:10]:
            rows.append(
                {"gene": g, "cytokine": "CYT_PARTIAL", "celltype": "A_ref", "log_fc": 2.0, "adj_p_value": 0.001}
            )
        # Add genes not in adata
        for i in range(5):
            rows.append(
                {
                    "gene": f"MISSING{i}",
                    "cytokine": "CYT_PARTIAL",
                    "celltype": "A_ref",
                    "log_fc": 2.0,
                    "adj_p_value": 0.001,
                }
            )
        df = pd.DataFrame(rows)

        gsd, gsdf = _get_genesets(
            adata, df, celltype_signature="A_ref", direction="upregulated", threshold_pval=0.01, threshold_lfc=1.0
        )
        assert len(gsd["CYT_PARTIAL"]) == 10
        row = gsdf[gsdf["cytokine"] == "CYT_PARTIAL"].iloc[0]
        assert row["num_genes_signature"] == 15
        assert row["num_shared_genes_signature"] == 10
        assert abs(row["frac_shared_genes_signature"] - 10 / 15) < 1e-6


# ---------------------------------------------------------------------------
# Tests: _compute_mu_and_sigma
# ---------------------------------------------------------------------------


class TestComputeMuAndSigma:
    def test_returns_expected_keys(self, adata):
        result = _compute_mu_and_sigma(adata, contrast_column="disease_state", condition="condition_1")
        assert set(result.keys()) == {"mu", "sigma", "num_cells"}

    def test_num_cells_correct(self, adata):
        result = _compute_mu_and_sigma(adata, contrast_column="disease_state", condition="condition_1")
        # 50 cells of type A + 50 cells of type B in condition_1
        assert result["num_cells"] == 100

    def test_mu_shape(self, adata, gene_names):
        result = _compute_mu_and_sigma(adata, contrast_column="disease_state", condition="condition_1")
        assert result["mu"].shape == (len(gene_names),)
        assert result["sigma"].shape == (len(gene_names),)

    def test_known_values(self):
        """Verify mu and sigma against manually computed values."""
        X = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], dtype=np.float32)
        obs = pd.DataFrame({"cond": ["A", "A", "A"]}, index=["c0", "c1", "c2"])
        adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=["g0", "g1"]))

        result = _compute_mu_and_sigma(adata, contrast_column="cond", condition="A")
        np.testing.assert_allclose(result["mu"], [3.0, 4.0], atol=1e-5)
        np.testing.assert_allclose(result["sigma"], [2.0, 2.0], atol=1e-5)
        assert result["num_cells"] == 3


# ---------------------------------------------------------------------------
# Tests: _compute_s2n
# ---------------------------------------------------------------------------


class TestComputeS2N:
    def test_output_shape(self, adata, gene_names):
        stats, num_cells = _compute_s2n(
            adata, contrast_column="disease_state", condition_1="condition_1", condition_2="condition_2"
        )
        assert stats.shape == (len(gene_names), 1)
        assert "condition_1_vs_condition_2" in stats.columns
        assert num_cells.shape == (1, 2)

    def test_known_s2n(self):
        """Manually verify S2N = (mu1 - mu2) / (sigma1 + sigma2 + eps)."""
        X = np.array(
            [
                [4.0, 1.0],
                [6.0, 3.0],
                [1.0, 4.0],
                [3.0, 6.0],
            ],
            dtype=np.float32,
        )
        obs = pd.DataFrame({"cond": ["A", "A", "B", "B"]}, index=[f"c{i}" for i in range(4)])
        adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=["g0", "g1"]))

        stats, num_cells = _compute_s2n(adata, contrast_column="cond", condition_1="A", condition_2="B")

        mu_a = np.array([5.0, 2.0])
        mu_b = np.array([2.0, 5.0])
        sigma_a = np.std([4.0, 6.0], ddof=1), np.std([1.0, 3.0], ddof=1)
        sigma_b = np.std([1.0, 3.0], ddof=1), np.std([4.0, 6.0], ddof=1)

        expected_s2n = (mu_a - mu_b) / (np.array(sigma_a) + np.array(sigma_b) + 1e-8)
        np.testing.assert_allclose(stats.values.flatten(), expected_s2n, atol=1e-5)
        assert num_cells.loc["A_vs_B", "num_cells_1"] == 2
        assert num_cells.loc["A_vs_B", "num_cells_2"] == 2

    def test_precomputed_stats_matches_direct(self, adata):
        """Passing precomputed stats should give identical results."""
        stats_direct, _ = _compute_s2n(
            adata, contrast_column="disease_state", condition_1="condition_1", condition_2="condition_2"
        )
        precomputed = {
            "condition_1": _compute_mu_and_sigma(adata, "disease_state", "condition_1"),
            "condition_2": _compute_mu_and_sigma(adata, "disease_state", "condition_2"),
        }
        stats_pre, _ = _compute_s2n(
            adata,
            contrast_column="disease_state",
            condition_1="condition_1",
            condition_2="condition_2",
            precomputed_stats=precomputed,
        )
        pd.testing.assert_frame_equal(stats_direct, stats_pre)

    def test_boosted_genes_positive_s2n(self, adata):
        """First 30 genes are boosted in condition_1 for cell type A.

        When using ALL cells, the effect is diluted by cell type B, but
        should still be positive for the boosted genes.
        """
        stats, _ = _compute_s2n(
            adata, contrast_column="disease_state", condition_1="condition_1", condition_2="condition_2"
        )
        # First 30 genes should have positive S2N (condition_1 higher than condition_2)
        assert (stats.iloc[:30].values > 0).all()


# ---------------------------------------------------------------------------
# Tests: _compute_ranking_statistic
# ---------------------------------------------------------------------------


class TestComputeRankingStatistic:
    def test_single_contrast(self, adata):
        contrasts = [("condition_1", "condition_2")]
        rnk_stats, num_cells = _compute_ranking_statistic(adata, "disease_state", contrasts)
        assert rnk_stats.shape[1] == 1
        assert num_cells.shape[0] == 1

    def test_multiple_contrasts(self, adata):
        contrasts = [("condition_1", "condition_2"), ("condition_2", "condition_1")]
        rnk_stats, num_cells = _compute_ranking_statistic(adata, "disease_state", contrasts)
        assert rnk_stats.shape[1] == 2
        assert num_cells.shape[0] == 2

    def test_reverse_contrast_is_negated(self, adata):
        contrasts = [("condition_1", "condition_2"), ("condition_2", "condition_1")]
        rnk_stats, _ = _compute_ranking_statistic(adata, "disease_state", contrasts)
        np.testing.assert_allclose(
            rnk_stats.iloc[:, 0].values,
            -rnk_stats.iloc[:, 1].values,
            atol=1e-5,
        )


# ---------------------------------------------------------------------------
# Tests: run_one_enrichment_test (integration)
# ---------------------------------------------------------------------------


class TestRunOneEnrichmentTest:
    def test_basic_hcd_run(self, adata, hcd_df):
        result = run_one_enrichment_test(
            adata=adata,
            df=hcd_df,
            celltype_combo=("A", "A_ref"),
            celltype_column="cell_type",
            contrasts_combo=("condition_1", "condition_2"),
            contrast_column="disease_state",
            direction="upregulated",
            threshold_lfc=1.0,
            threshold_expression=0.0,
            threshold_pval=0.01,
            min_size=5,
            permutation_num=100,
            seed=42,
            threads=1,
        )
        assert result is not None
        assert isinstance(result, pd.DataFrame)
        assert "cytokine" in result.columns
        assert "NES" in result.columns
        assert "FDR q-val" in result.columns

    def test_hcd_result_columns(self, adata, hcd_df):
        result = run_one_enrichment_test(
            adata=adata,
            df=hcd_df,
            celltype_combo=("A", "A_ref"),
            celltype_column="cell_type",
            contrasts_combo=("condition_1", "condition_2"),
            contrast_column="disease_state",
            direction="upregulated",
            threshold_lfc=1.0,
            threshold_expression=0.0,
            threshold_pval=0.01,
            min_size=5,
            permutation_num=100,
            seed=42,
            threads=1,
        )
        expected_cols = {
            "cytokine",
            "NES",
            "FDR q-val",
            "contrast",
            "celltype_adata",
            "celltype_signature",
            "celltype_combo",
            "direction",
            "threshold_pval",
            "threshold_lfc",
            "threshold_expression",
            "num_cells_1",
            "num_cells_2",
            "num_genes_signature",
            "num_shared_genes_signature",
            "frac_shared_genes_signature",
        }
        assert expected_cols.issubset(result.columns)

    def test_hcd_metadata_values(self, adata, hcd_df):
        result = run_one_enrichment_test(
            adata=adata,
            df=hcd_df,
            celltype_combo=("A", "A_ref"),
            celltype_column="cell_type",
            contrasts_combo=("condition_1", "condition_2"),
            contrast_column="disease_state",
            direction="upregulated",
            threshold_lfc=1.0,
            threshold_expression=0.0,
            threshold_pval=0.01,
            min_size=5,
            permutation_num=100,
            seed=42,
            threads=1,
        )
        assert (result["celltype_adata"] == "A").all()
        assert (result["celltype_signature"] == "A_ref").all()
        assert (result["celltype_combo"] == "A (A_ref)").all()
        assert (result["direction"] == "upregulated").all()
        assert (result["contrast"] == "condition_1_vs_condition_2").all()

    def test_missing_celltype_returns_none(self, adata, hcd_df):
        result = run_one_enrichment_test(
            adata=adata,
            df=hcd_df,
            celltype_combo=("NONEXISTENT", "A_ref"),
            celltype_column="cell_type",
            contrasts_combo=("condition_1", "condition_2"),
            contrast_column="disease_state",
            direction="upregulated",
            min_size=5,
            permutation_num=100,
            seed=42,
            threads=1,
        )
        assert result is None

    def test_cip_run(self, adata, cip_df):
        result = run_one_enrichment_test(
            adata=adata,
            df=cip_df,
            celltype_combo=("A", "A_ref"),
            celltype_column="cell_type",
            contrasts_combo=("condition_1", "condition_2"),
            contrast_column="disease_state",
            direction="upregulated",
            min_size=5,
            permutation_num=100,
            seed=42,
            threads=1,
        )
        assert result is not None
        assert "CIP" in result.columns
        assert (result["direction"] == "upregulated").all()

    def test_custom_run(self, adata, custom_df):
        result = run_one_enrichment_test(
            adata=adata,
            df=custom_df,
            celltype_combo=("A", "A_ref"),
            celltype_column="cell_type",
            contrasts_combo=("condition_1", "condition_2"),
            contrast_column="disease_state",
            direction="upregulated",
            min_size=5,
            permutation_num=100,
            seed=42,
            threads=1,
        )
        assert result is not None
        assert "query_program" in result.columns
        assert (result["direction"] == "custom input").all()

    def test_multiple_contrasts(self, adata, hcd_df):
        result = run_one_enrichment_test(
            adata=adata,
            df=hcd_df,
            celltype_combo=("A", "A_ref"),
            celltype_column="cell_type",
            contrasts_combo=[("condition_1", "condition_2"), ("condition_2", "condition_1")],
            contrast_column="disease_state",
            direction="upregulated",
            threshold_lfc=1.0,
            threshold_expression=0.0,
            threshold_pval=0.01,
            min_size=5,
            permutation_num=100,
            seed=42,
            threads=1,
        )
        assert result is not None
        assert set(result["contrast"].unique()) == {"condition_1_vs_condition_2", "condition_2_vs_condition_1"}

    def test_reproducibility_with_seed(self, adata, hcd_df):
        kwargs = {
            "adata": adata,
            "df": hcd_df,
            "celltype_combo": ("A", "A_ref"),
            "celltype_column": "cell_type",
            "contrasts_combo": ("condition_1", "condition_2"),
            "contrast_column": "disease_state",
            "direction": "upregulated",
            "threshold_lfc": 1.0,
            "threshold_expression": 0.0,
            "threshold_pval": 0.01,
            "min_size": 5,
            "permutation_num": 100,
            "seed": 42,
            "threads": 1,
        }
        result1 = run_one_enrichment_test(**kwargs)
        result2 = run_one_enrichment_test(**kwargs)
        pd.testing.assert_frame_equal(result1, result2)

    def test_boosted_genes_show_enrichment(self, adata, hcd_df):
        """CYT1 (first 20 genes) overlaps heavily with the boosted region.

        In cell type A, condition_1 has higher expression for those genes,
        so CYT1 should show positive NES for condition_1 vs condition_2.
        """
        result = run_one_enrichment_test(
            adata=adata,
            df=hcd_df,
            celltype_combo=("A", "A_ref"),
            celltype_column="cell_type",
            contrasts_combo=("condition_1", "condition_2"),
            contrast_column="disease_state",
            direction="upregulated",
            threshold_lfc=1.0,
            threshold_expression=0.0,
            threshold_pval=0.01,
            min_size=5,
            permutation_num=100,
            seed=42,
            threads=1,
        )
        cyt1_nes = result.loc[result["cytokine"] == "CYT1", "NES"].values[0]
        assert float(cyt1_nes) > 0, f"Expected positive NES for CYT1, got {cyt1_nes}"


# ---------------------------------------------------------------------------
# Tests: run_all_enrichment_test (integration)
# ---------------------------------------------------------------------------


class TestRunAllEnrichmentTest:
    def test_basic_run(self, adata, hcd_df):
        result = run_all_enrichment_test(
            adata=adata,
            df=hcd_df,
            celltype_combos=[("A", "A_ref")],
            celltype_column="cell_type",
            contrasts_combo=("condition_1", "condition_2"),
            contrast_column="disease_state",
            direction="upregulated",
            threshold_lfc=1.0,
            threshold_expression=0.0,
            threshold_pval=0.01,
            min_size=5,
            permutation_num=100,
            seed=42,
            threads=1,
        )
        assert result is not None
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0

    def test_multiple_thresholds_expand_results(self, adata, hcd_df):
        result = run_all_enrichment_test(
            adata=adata,
            df=hcd_df,
            celltype_combos=[("A", "A_ref")],
            celltype_column="cell_type",
            contrasts_combo=("condition_1", "condition_2"),
            contrast_column="disease_state",
            direction="upregulated",
            threshold_lfc=[0.5, 1.0],
            threshold_expression=[0.0, 0.1],
            threshold_pval=0.01,
            min_size=5,
            permutation_num=100,
            seed=42,
            threads=1,
        )
        # With 2 lfc thresholds x 2 expression thresholds = 4 runs
        # Each run produces results for each cytokine that passes min_size
        assert len(result["threshold_lfc"].unique()) == 2
        assert len(result["threshold_expression"].unique()) == 2

    def test_multiple_celltypes(self, adata, hcd_df, gene_names):
        # Add cell type B entries to hcd
        rows_b = [
            {"gene": g, "cytokine": "CYT_B", "celltype": "B_ref", "log_fc": 2.0, "adj_p_value": 0.001}
            for g in gene_names[:15]
        ]
        hcd_extended = pd.concat([hcd_df, pd.DataFrame(rows_b)], ignore_index=True)

        result = run_all_enrichment_test(
            adata=adata,
            df=hcd_extended,
            celltype_combos=[("A", "A_ref"), ("B", "B_ref")],
            celltype_column="cell_type",
            contrasts_combo=("condition_1", "condition_2"),
            contrast_column="disease_state",
            direction="upregulated",
            threshold_lfc=1.0,
            threshold_expression=0.0,
            threshold_pval=0.01,
            min_size=5,
            permutation_num=100,
            seed=42,
            threads=1,
        )
        assert set(result["celltype_adata"].unique()) == {"A", "B"}
