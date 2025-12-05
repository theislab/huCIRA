from anndata import AnnData 
import pandas as pd
import scanpy as sc

def load_cytokine_dict_data():
  """
  To be changed: Currently referring to local path.
  """
  cytokine_dict = pd.read_csv("/home/icb/jenni.liu/projects/cytokine_dict_folder/DEGs_all_with_pbs_DE_count.csv", index_col=0)
  return cytokine_dict


def load_MS_data():
  """
  To be changed. Currently referring to local path.
  """
  adata = sc.read_h5ad("/home/icb/jenni.liu/projects/cytokine_dict_folder/Schafflick20_MS_CSF.h5ad")
  return adata


def load_cytokine_info():
    """
    Cytokine information sheet includes information about sender and receptor genes (for cell-cell communication plot).
    """

    url = "https://github.com/theislab/cytokine_dict/edit/main/src/cytokine_dict/data/20250125_cytokine_info_with_functional_classification_LV.xlsx"
    cytokine_info = pd.read_excel(url, sheet_name="all_cytokines")
    return cytokine_info
