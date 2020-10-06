import pandas as pd 
import re
from maayanlab_bioinformatics.harmonization import ncbi_genes


def get_ensembl_id(ids):
    ids = "".join(ids)
    ensembl = re.findall("Ensembl:([A-Z]*[0-9]{11})", ids)
    if len(ensembl) > 0: return ensembl[0]
    return None 
       

def get_symbol(df): 
    """
    Convert Ensembl ID to gene symbol for the index (column) of DataFrame.
    """
    ncbi = pd.DataFrame(ncbi_genes.ncbi_genes_fetch())
    dbxrefs = ncbi.loc[:, 'dbXrefs']
    symbols = ncbi.loc[:, 'Symbol']
    eid_to_symbol = {}
    for i in range(len(dbxrefs)): 
        item = dbxrefs[i]
        symbol = symbols[i]
        ensembl = get_ensembl_id(item)
        if ensembl: eid_to_symbol[ensembl] = symbol
    def id_to_symbol(key):
        if (key in eid_to_symbol):
            return eid_to_symbol[key]
        else:
            return key
    data_ensembl_ids = df.index.to_list()
    data_symbols = [ id_to_symbol(key) for key in data_ensembl_ids ]
    df.index = data_symbols
    return df