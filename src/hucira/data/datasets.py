import pandas as pd
import scanpy as sc
import os
import requests



def load_human_cytokine_dict():
    """To be changed: Currently referring to local path."""
    cytokine_dict = pd.read_csv(
        "/home/icb/jenni.liu/projects/cytokine_dict_folder/DEGs_all_with_pbs_DE_count.csv", index_col=0
    )
    return cytokine_dict


def load_MS_CSF_data(save_dir="",
                     force_download=False):
    """
    Download and load the MS dataset from automatically.
    Xu, Chenling (2021). MS_CSF.h5ad. figshare. Dataset. https://doi.org/10.6084/m9.figshare.14356661.v1
    
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
    
    url = "https://figshare.com/ndownloader/files/27405182"
    if save_dir == "":
        save_dir = os.getcwd()

    local_path = os.path.join(save_dir, "MS_CSF.h5ad")

    # Download only if not already in cache
    if force_download or not os.path.exists(local_path):
        print("Downloading lupus dataset from CELLxGENE...")
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(local_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        print(f"Download complete: {local_path}")
    else:
        print(f"Using cached file: {local_path}")

    # Load with scanpy
    return sc.read_h5ad(local_path)



def load_Lupus_data(save_dir="",
                    force_download=False):
    """
    Download and load the lupus dataset from CELLxGENE automatically.
    Richard K. Perez et al. ,Single-cell RNA-seq reveals cell type–specific molecular and genetic associations to lupus.Science376,eabf1970(2022).DOI:10.1126/science.abf1970

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
    local_path = os.path.join(save_dir, "lupus.h5ad")

    # Download only if not already in cache
    if force_download or not os.path.exists(local_path):
        print("Downloading lupus dataset from CELLxGENE...")
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(local_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        print(f"Download complete: {local_path}")
    else:
        print(f"Using cached file: {local_path}")

    # Load with scanpy
    return sc.read_h5ad(local_path)


def load_cytokine_info():
    """Cytokine information sheet includes information about sender and receptor genes (for cell-cell communication plot)."""
    url = "https://github.com/theislab/cytokine_dict/edit/main/src/cytokine_dict/data/20250125_cytokine_info_with_functional_classification_LV.xlsx"
    cytokine_info = pd.read_excel(url, sheet_name="all_cytokines")
    return cytokine_info
