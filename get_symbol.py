import pandas as pd 
import re
from maayanlab_bioinformatics.harmonization import ncbi_genes


def get_ensembl_id(ids):
    ids = "".join(ids)
    ensembl = re.findall("Ensembl:(.*)", ids)
    if (len(ensembl) == 1):
        return ensembl[0]
    else:
        return None
       

def get_symbol(df): 
    """
    Convert Ensembl ID to gene symbol for the index (column) of DataFrame.
    """
    ncbi = pd.DataFrame(ncbi_genes.ncbi_genes_fetch())
    all_ids = ncbi.dbXrefs.values
    ensembl_ids = [ get_ensembl_id(ids) for ids in all_ids]
    ncbi = ncbi[["dbXrefs", "Symbol", "type_of_gene"]]
    ncbi["ensembl"] = ensembl_ids
    ncbi = ncbi.drop(columns=["dbXrefs"])
    ncbi = ncbi.set_index("ensembl")
    def id_to_symbol(key):
        if (key in ensembl_to_symbol):
            return ensembl_to_symbol[key]
        else:
            return key
    ensembl_to_symbol = ncbi.to_dict()["Symbol"]
    data_ensembl_ids = df.index.to_list()
    data_symbols = [ id_to_symbol(key) for key in data_ensembl_ids ]
    df.index = data_symbols
    df.columns = data_symbols
    return df