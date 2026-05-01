import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from hucira.tl.enrichment_test import _compute_log_fold_change, run_one_enrichment_test


@pytest.fixture
def synthetic_data():
    """Create a synthetic dataset for enrichment testing.

    - 200 genes, 100 cells
    - Two cell types: "T cell" (query) mapped to "T" (signature)
    - Two conditions: "disease" vs "healthy"
    - A custom gene program with 15 genes so it passes min_size=10
    """
    rng = np.random.default_rng(42)
    n_cells = 100
    n_genes = 200
    gene_names = [f"gene_{i}" for i in range(n_genes)]

    # Create expression matrix with condition-dependent signal in first 15 genes
    X = rng.poisson(lam=2, size=(n_cells, n_genes)).astype(np.float32)
    # Boost first 15 genes in "disease" cells (first 50 cells)
    X[:50, :15] += 3.0

    adata = ad.AnnData(X=X)
    adata.var_names = gene_names
    adata.obs["cell_type"] = "T cell"
    adata.obs["condition"] = ["disease"] * 50 + ["healthy"] * 50

    # Gene program DataFrame (query_program format)
    program_genes = gene_names[:15]
    df = pd.DataFrame({"gene": program_genes, "query_program": "program_A", "celltype": "T"})

    return adata, df


def test_custom_ranking_fn_runs(synthetic_data):
    """run_one_enrichment_test completes with a custom ranking function."""
    adata, df = synthetic_data

    result = run_one_enrichment_test(
        adata=adata,
        df=df,
        celltype_combo=("T cell", "T"),
        contrast_combo=("disease", "healthy"),
        celltype_column="cell_type",
        contrast_column="condition",
        direction="upregulated",
        threshold_lfc=0.0,
        threshold_pval=1.0,
        min_size=5,
        permutation_num=100,
        threads=1,
        ranking_statistic_fn=_compute_log_fold_change,
    )

    assert result is not None
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0
    assert "NES" in result.columns
    assert "contrast" in result.columns


def test_default_ranking_fn_runs(synthetic_data):
    """run_one_enrichment_test completes with default S2N (ranking_statistic_fn=None)."""
    adata, df = synthetic_data

    result = run_one_enrichment_test(
        adata=adata,
        df=df,
        celltype_combo=("T cell", "T"),
        contrast_combo=("disease", "healthy"),
        celltype_column="cell_type",
        contrast_column="condition",
        direction="upregulated",
        threshold_lfc=0.0,
        threshold_pval=1.0,
        min_size=5,
        permutation_num=100,
        threads=1,
        ranking_statistic_fn=None,
    )

    assert result is not None
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0


def test_custom_ranking_fn_with_sparse_matrix(synthetic_data):
    """Custom ranking function works with sparse matrices."""
    adata, df = synthetic_data
    adata.X = sp.csr_matrix(adata.X)

    result = run_one_enrichment_test(
        adata=adata,
        df=df,
        celltype_combo=("T cell", "T"),
        contrast_combo=("disease", "healthy"),
        celltype_column="cell_type",
        contrast_column="condition",
        direction="upregulated",
        threshold_lfc=0.0,
        threshold_pval=1.0,
        min_size=5,
        permutation_num=100,
        threads=1,
        ranking_statistic_fn=_compute_log_fold_change,
    )

    assert result is not None
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0
