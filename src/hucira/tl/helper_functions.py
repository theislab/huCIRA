import warnings


def create_celltype_combos(
    adata_celltypes: list[str] | None = None, hcd_celltypes: list[str] | None = None
) -> list[tuple[str, str]] | None:
    """Reformat cell-type names for :func:`run_one_enrichment_test`.

    Creates a ``celltype_combos`` list of tuples that matches cell types of the
    query *adata* to the closest available cell types from the human cytokine
    dictionary.  The lists and their order have to be manually prepared to
    fit/match.  This is just a helper -- the list can be manually curated.

    Parameters
    ----------
    adata_celltypes : list of str or None
        Ordered cell types of interest from the adata object (matching
        *hcd_celltypes*).
    hcd_celltypes : list of str or None
        Ordered cell types of interest from the human cytokine dictionary
        (matching *adata_celltypes*).

    Returns
    -------
    list of tuple
        A list of matched cell-type pairs, e.g.
        ``[('B cell', 'B'), ('CD4-positive, alpha-beta T cell', 'CD4'), ...]``.
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

    celltype_combos = list(zip(adata_celltypes, hcd_celltypes, strict=True))
    return celltype_combos
