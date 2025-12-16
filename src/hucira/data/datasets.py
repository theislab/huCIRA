import pandas as pd
import scanpy as sc
import os
import requests



def load_human_cytokine_dict(save_dir="", force_download=False):
    """
    Download and load our Human Cytokine Dictionary from Parse Biosciences:
    https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/
    
    Parameters
    ----------
    save_dir : str
        Directory where the file will be saved.
    force_download : bool
        Allows user to force a fresh download 

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
        print("Downloading Human Cytokine Dictionary from Parse Biosciences...")
        cytokine_dict = pd.read_csv(url, index_col=0)
        cytokine_dict.to_csv(local_path)
    else:
        print(f"Loading from: {local_path}")
        cytokine_dict = pd.read_csv(local_path, index_col=0)

    return cytokine_dict


def load_MS_CSF_data(save_dir="", force_download=False):
    """
    Download and load the MS dataset automatically.
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
        MS adata object.
    """
    
    url = "https://figshare.com/ndownloader/files/27405182"
    if save_dir == "":
        save_dir = os.getcwd()     
    os.makedirs(save_dir, exist_ok=True)
    local_path = os.path.join(save_dir, "MS_CSF.h5ad")

    # Download only if not already in directory
    if force_download or not os.path.exists(local_path):
        print("Downloading MS dataset from figshare...")
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(local_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        print(f"Download complete: {local_path}")
    else:
        print(f"Loading from: {local_path}")

    # Load with scanpy
    return sc.read_h5ad(local_path)



def load_Lupus_data(save_dir="", force_download=False):
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
    os.makedirs(save_dir, exist_ok=True)
    local_path = os.path.join(save_dir, "lupus.h5ad")

    # Download only if not already in directory
    if force_download or not os.path.exists(local_path):
        print("Downloading lupus dataset from CELLxGENE...")
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(local_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        print(f"Download complete: {local_path}")
    else:
        print(f"Loading from: {local_path}")

    # Load with scanpy
    return sc.read_h5ad(local_path)


def load_cytokine_info(save_dir="", force_download=False):
    """
    Download and load Cytokine information sheet: includes information about sender and receptor genes (for cell-cell communication plot)
    
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
        print("Downloading Cytokine Information sheet...")
        cytokine_info = pd.read_excel(url, sheet_name="all_cytokines", engine='openpyxl')
        cytokine_info.to_excel(local_path, sheet_name="all_cytokines")
    else:
        print(f"Loading from: {local_path}")
        cytokine_info = pd.read_excel(local_path, index_col=0)

    return cytokine_info




def load_CIP_signatures(save_dir="", force_download=False):
    """
    Download and load metadata file (sheet "13.CIP_activations") from supplemental data: information about CIPs (cytokine induced gene programs).
    
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

    url = "https://www.biorxiv.org/content/biorxiv/early/2025/12/15/2025.12.12.693897/DC2/embed/media-2.xlsx?download=true"
    if save_dir == "":
        save_dir = os.getcwd()
    os.makedirs(save_dir, exist_ok=True)
    local_path = os.path.join(save_dir, "CIP_signatures.xlsx")

    if force_download or not os.path.exists(local_path):
        print("Downloading Cytokine induced gene programs sheet...")
        CIP_signatures = pd.read_excel(url, sheet_name="13.CIP_activations", engine='openpyxl')
        CIP_signatures.to_excel(local_path, sheet_name="13.CIP_activations", index=False)
    else:
        print(f"Loading from: {local_path}")
        CIP_signatures = pd.read_excel(local_path, index_col=0)

    return CIP_signatures

    
