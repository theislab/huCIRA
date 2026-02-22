import contextlib
import os

import pandas as pd
import requests
import requests_cache
import scanpy as sc
from anndata import AnnData
from tqdm.auto import tqdm


def _download_file(url: str, local_path: str, description: str, disable_cache: bool = False) -> None:
    """Download a file with a tqdm progress bar."""
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


def load_human_cytokine_dict(save_dir="", force_download=False, exclude_well_biased_genes=True) -> pd.DataFrame:
    """
    Download and load our Human Cytokine Dictionary from Parse Biosciences:

    Source: https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/

    Parameters
    ----------
    save_dir : str
        Directory where the file will be saved.
    force_download : bool
        Allows user to force a fresh download
    exclude_well_biased_genes : bool
        If True, exclude genes that are well biased according to out analysis in the original publication.


    Returns
    -------
    cytokine_dict : pandas.DataFrame
        Human Cytokine Dictionary adata object.
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
        print(f"Loading from: {local_path}")
        cytokine_dict = pd.read_csv(local_path, index_col=0)
        cytokine_dict = cytokine_dict.reset_index(drop=True)

    if exclude_well_biased_genes:
        cytokine_dict = cytokine_dict.loc[~cytokine_dict.well_biased]

    # Make sure that these three cytokines are removed from dictionary if object not updated yet online.
    revision_cytokines = ["TGF-beta1", "IL-18", "C3a"]
    cytokine_dict = cytokine_dict[~cytokine_dict["cytokine"].isin(revision_cytokines)]
    cytokine_dict = cytokine_dict.reset_index(drop=True)

    return cytokine_dict


def load_MS_CSF_data(save_dir="", force_download=False) -> AnnData:
    """
    Download and load the MS dataset automatically.

    Reference: Xu, Chenling (2021). MS_CSF.h5ad. figshare. Dataset. https://doi.org/10.6084/m9.figshare.14356661.v1

    Parameters
    ----------
    save_dir : str
        Directory where the file will be saved.
    force_download : bool
        Allows user to force a fresh download from CellxGene

    Returns
    -------
    adata : AnnData
        MS adata object.
    """
    url = "https://ndownloader.figshare.com/files/27405182"

    if save_dir == "":
        save_dir = os.getcwd()
    os.makedirs(save_dir, exist_ok=True)

    local_path = os.path.join(save_dir, "MS_CSF.h5ad")

    # Download only if not already in directory
    if force_download or not os.path.exists(local_path):
        _download_file(url, local_path, "MS dataset (figshare)", disable_cache=True)
    else:
        print(f"Loading from: {local_path}")

    # Load with scanpy
    return sc.read_h5ad(local_path)


def load_Lupus_data(save_dir="", force_download=False) -> AnnData:
    """
    Download and load the lupus dataset from CELLxGENE automatically.

    Reference: Richard K. Perez et al. ,Single-cell RNA-seq reveals cell type–specific molecular and genetic associations to lupus.Science376,eabf1970(2022).DOI:10.1126/science.abf1970

    Parameters
    ----------
    save_dir : str
        Directory where the file will be saved.
    force_download : bool
        Allows user to force a fresh download from CellxGene

    Returns
    -------
    adata : AnnData
        Lupus adata object.
    """
    url = "https://datasets.cellxgene.cziscience.com/4118e166-34f5-4c1f-9eed-c64b90a3dace.h5ad"
    if save_dir == "":
        save_dir = os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    local_path = os.path.join(save_dir, "lupus.h5ad")

    # Download only if not already in directory
    if force_download or not os.path.exists(local_path):
        _download_file(url, local_path, "Lupus dataset (CELLxGENE)")
    else:
        print(f"Loading from: {local_path}")

    # Load with scanpy
    return sc.read_h5ad(local_path)


def load_cytokine_info(save_dir="", force_download=False) -> pd.DataFrame:
    """
    Download and load Cytokine information sheet.

    The information sheet includes information about sender and receptor genes (for cell-cell communication plot).

    Parameters
    ----------
    save_dir : str
        Directory where the file will be saved.
    force_download : bool
        Allows user to force a fresh download

    Returns
    -------
    cytokine_info : pandas.DataFrame
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
        print(f"Loading from: {local_path}")
        cytokine_info = pd.read_excel(local_path)

    return cytokine_info


def load_CIP_signatures(save_dir="", force_download=False) -> pd.DataFrame:
    """
    Download and load metadata file (sheet "13.CIP_activations") from supplemental data.

    This file contains information about CIPs (cytokine induced gene programs).

    Parameters
    ----------
    save_dir : str
        Directory where the file will be saved.
    force_download : bool
        Allows user to force a fresh download

    Returns
    -------
    CIP_signatures : pandas.DataFrame
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
        print(f"Loading from: {local_path}")
        CIP_signatures = pd.read_csv(local_path, index_col=0)

    return CIP_signatures
