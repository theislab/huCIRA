import warnings


def create_celltype_combos(adata_celltypes=None, hcd_celltypes=None):
    """Reformatting celltype of input data for run_enrichment_test()

    Creates celltype_combos tuple, which is container that matches celltypes of query adata to the closest fit available celltypes from the human cytokine dictionary. The lists and their order have to be manually prepared to fit/match.
    This is just a helper, the tuple can be manually curated.

    Parameters
    ----------
    adata_celltypes
        ordered celltypes of interest from adata object (matching with hcd_celltypes)
    hcd_celltypes
        ordered celltypes of interest from human cytokine dictionary (matching with adata_celltypes)

    Returns
    -------
    celltype_combo
        A tuple of matched celltype sets. Example output:
        celltype_combo = (('B cell', 'B'),
                         ('CD4-positive, alpha-beta T cell', 'CD4'),
                         ('CD8-positive, alpha-beta T cell', 'CD8'),
                         ('classical monocyte', 'CD14_Mono'),
                         ('non-classical monocyte', 'CD16_Mono'),
                         ('conventional dendritic cell', 'cDC'),
                         ('natural killer cell', 'NK_CD56hi'),
                         ('natural killer cell', 'NK_CD56low'))
    """
    # Make sure they are lists
    adata_celltypes = list(adata_celltypes) if adata_celltypes is not None else []
    hcd_celltypes = list(hcd_celltypes) if hcd_celltypes is not None else []

    if len(adata_celltypes) != len(hcd_celltypes):
        warnings.warn(
            f"The lists have different lengths: adata_celltypes={len(adata_celltypes)}, "
            f"hcd_celltypes={len(hcd_celltypes)}.\nCheck input lists. Did not return celltype_combos.",
            stacklevel=2,  # points to the caller of the function
        )
        return

    celltype_combos = tuple(zip(adata_celltypes, hcd_celltypes, strict=True))
    return celltype_combos
