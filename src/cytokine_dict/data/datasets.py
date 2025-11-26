from anndata import AnnData
import pandas as pd

def load_cytokine_dict_data():
  """
  To be changed. 
  """
  cytokine_dict = pd.read_csv("/home/icb/jenni.liu/projects/cytokine_dict_folder/DEGs_all_with_pbs_DE_count.csv", index_col=0))
  return cytokine_dict


def load_MS_data():
  """
  To be changed. 
  """
  adata = sc.read_h5ad("/home/icb/jenni.liu/projects/cytokine_dict_folder/Schafflick20_MS_CSF.h5ad")
  return adata
