import contextlib
import logging
import os

import pandas as pd
import requests
import requests_cache
import scanpy as sc
from anndata import AnnData
from tqdm.auto import tqdm

logger = logging.getLogger(__name__)


def _download_file(url: str, local_path: str, description: str, disable_cache: bool = False) -> None:
    """Download a file with a tqdm progress bar."""
    logger.info("Downloading %s from %s", description, url)
    ctx = requests_cache.disabled() if disable_cache else contextlib.nullcontext()
    with ctx:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            total_size = int(r.headers.get("content-length", 0))
            with tqdm(total=total_size, unit="B", unit_scale=True, desc=description) as pbar:
                with open(local_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
                        pbar.update(len(chunk))


def load_human_cytokine_dict(
    save_dir: str = "", force_download: bool = False, exclude_well_biased_genes: bool = True
) -> pd.DataFrame:
    """Download and load the Human Cytokine Dictionary.

    Source: https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/

    Parameters
    ----------
    save_dir : str
        Directory where the file will be saved.
    force_download : bool
        Allows user to force a fresh download.
    exclude_well_biased_genes : bool
        If True, exclude genes that are well biased according to our
        analysis in the original publication.

    Returns
    -------
    cytokine_dict : pandas.DataFrame
        Human Cytokine Dictionary as a DataFrame.
    """
    url = "https://cdn.parsebiosciences.com/gigalab/10m/DEGs.csv"
    if save_dir == "":
        save_dir = os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    local_path = os.path.join(save_dir, "human_cytokine_dict.csv")

    if force_download or not os.path.exists(local_path):
        _download_file(url, local_path, "Human Cytokine Dictionary")
        cytokine_dict = pd.read_csv(local_path, index_col=0)
        cytokine_dict = cytokine_dict.reset_index(drop=True)
        cytokine_dict.to_csv(local_path)
    else:
        logger.info("Loading from: %s", local_path)
        cytokine_dict = pd.read_csv(local_path, index_col=0)
        cytokine_dict = cytokine_dict.reset_index(drop=True)

    if exclude_well_biased_genes:
        cytokine_dict = cytokine_dict.loc[~cytokine_dict.well_biased]

    # Make sure that these three cytokines are removed from dictionary if object not updated yet online.
    revision_cytokines = ["TGF-beta1", "IL-18", "C3a"]
    cytokine_dict = cytokine_dict[~cytokine_dict["cytokine"].isin(revision_cytokines)]
    cytokine_dict = cytokine_dict.reset_index(drop=True)

    return cytokine_dict


def load_multiple_sclerosis_data(save_dir: str = "", force_download: bool = False) -> AnnData:
    """Download and load the multiple sclerosis dataset.

    Reference: Xu, Chenling (2021). MS_CSF.h5ad. figshare. Dataset.
    https://doi.org/10.6084/m9.figshare.14356661.v1

    Parameters
    ----------
    save_dir : str
        Directory where the file will be saved.
    force_download : bool
        Allows user to force a fresh download from figshare.

    Returns
    -------
    adata : AnnData
        Multiple sclerosis AnnData object.
    """
    url = "https://ndownloader.figshare.com/files/27405182"

    if save_dir == "":
        save_dir = os.getcwd()
    os.makedirs(save_dir, exist_ok=True)

    local_path = os.path.join(save_dir, "MS_CSF.h5ad")

    # Download only if not already in directory
    if force_download or not os.path.exists(local_path):
        _download_file(url, local_path, "MS dataset", disable_cache=True)
    else:
        logger.info("Loading from: %s", local_path)

    # Load with scanpy
    return sc.read_h5ad(local_path)


def load_lupus_data(save_dir: str = "", force_download: bool = False) -> AnnData:
    """Download and load the lupus dataset from CELLxGENE.

    Reference: Perez et al., Single-cell RNA-seq reveals cell
    type-specific molecular and genetic associations to lupus.
    Science 376, eabf1970 (2022). DOI:10.1126/science.abf1970

    Parameters
    ----------
    save_dir : str
        Directory where the file will be saved.
    force_download : bool
        Allows user to force a fresh download from CellxGene.

    Returns
    -------
    adata : AnnData
        Lupus AnnData object.
    """
    url = "https://datasets.cellxgene.cziscience.com/4118e166-34f5-4c1f-9eed-c64b90a3dace.h5ad"
    if save_dir == "":
        save_dir = os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    local_path = os.path.join(save_dir, "lupus.h5ad")

    # Download only if not already in directory
    if force_download or not os.path.exists(local_path):
        _download_file(url, local_path, "Lupus dataset")
    else:
        logger.info("Loading from: %s", local_path)

    # Load with scanpy
    return sc.read_h5ad(local_path)


def load_cytokine_info(save_dir: str = "", force_download: bool = False) -> pd.DataFrame:
    """Download and load the cytokine information sheet.

    The information sheet includes sender and receptor genes used for the
    cell-cell communication plot.

    Parameters
    ----------
    save_dir : str
        Directory where the file will be saved.
    force_download : bool
        Allows user to force a fresh download.

    Returns
    -------
    cytokine_info : pandas.DataFrame
        Cytokine information including sender and receptor genes.
    """
    url = (
        "https://raw.githubusercontent.com/theislab/huCIRA/"
        "main/src/hucira/data/"
        "20250125_cytokine_info_with_functional_classification_LV.xlsx"
    )

    if save_dir == "":
        save_dir = os.getcwd()
    os.makedirs(save_dir, exist_ok=True)

    local_path = os.path.join(save_dir, "cytokine_info.xlsx")

    if force_download or not os.path.exists(local_path):
        _download_file(url, local_path, "Cytokine Information sheet")
        cytokine_info = pd.read_excel(local_path, sheet_name="all_cytokines", engine="openpyxl")
        cytokine_info.to_excel(local_path, sheet_name="all_cytokines")
    else:
        logger.info("Loading from: %s", local_path)
        cytokine_info = pd.read_excel(local_path)

    return cytokine_info


def load_CIP_signatures(save_dir: str = "", force_download: bool = False) -> pd.DataFrame:
    """Download and load CIP (cytokine-induced program) signatures.

    Parameters
    ----------
    save_dir : str
        Directory where the file will be saved.
    force_download : bool
        Allows user to force a fresh download.

    Returns
    -------
    CIP_signatures : pandas.DataFrame
        DataFrame containing CIP gene sets.
    """
    url = "https://raw.githubusercontent.com/theislab/huCIRA/main/src/hucira/data/df_cips_genesets.csv"
    if save_dir == "":
        save_dir = os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    local_path = os.path.join(save_dir, "CIP_signatures.csv")

    if force_download or not os.path.exists(local_path):
        _download_file(url, local_path, "CIP signatures")
        CIP_signatures = pd.read_csv(local_path, index_col=0)
        CIP_signatures.to_csv(local_path, index=False)
    else:
        logger.info("Loading from: %s", local_path)
        CIP_signatures = pd.read_csv(local_path, index_col=0)

    return CIP_signatures
