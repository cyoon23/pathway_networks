import numpy as np
import pandas as pd
import torch
from get_symbol import get_symbol
from prediction_pipeline import *
from preprocessing import preprocess


def load_data(expression_path, libraries, gr_truth_path): 
    norm = get_expression_data(expression_path)
    binary_matrix = get_binary_matrix(norm.index, libraries)
    gr_truth = get_gr_truth(norm.columns, gr_truth_path)
    list_inputs = prepare_data(norm, binary_matrix)
    return list_inputs, gr_truth

def get_expression_data(expression_path): 
    """ 
    Get expression data from the file path, turn Ensembl number
    indexing into gene symbols, normalize data. 
    """
    f = pd.read_csv(expression_path)
    f.index = f.iloc[:, 0] # Make ENSG genes as row indexing 
    f = f.iloc[:, 1:] # Remove first index column 
    # Normalize data 
    norm = preprocess(f)
    # Convert Ensembl number index to gene symbol
    norm = get_symbol(norm)
    return norm 

def get_binary_matrix(gene_expr, libraries):
    """ 
    Get binary matrix with genes as rows and pathways as columns. 
    If a gene is found in a given pathway, it is given a value of
    1. Else, 0. Only the list of genes in common between that found
    in the gene set libraries and the current dataset are used. 
    """
    function_to_genes = {}
    set_genes = set()
    for lib in libraries: 
        f2g, genes = gene_set_dictionaries(lib)
        function_to_genes.update(f2g)
        set_genes = set_genes | set(genes)
    common_genes = list(set_genes & set(gene_expr))
    binary_matrix = gs_binary_matrix(function_to_genes, set_genes).loc[common_genes]
    return binary_matrix

def get_gr_truth(list_samples, gr_truth_path):
    gr_truth = pd.read_csv(gr_truth_path, sep='\t')
    gr_truth = gr_truth.iloc[10, 1:]
    idx = [ "_".join(i.split("_")[:2]) for i in gr_truth.index ]
    gr_truth.index = idx
    gr_truth = pd.DataFrame(gr_truth).loc[list_samples]
    test = [ 1.0 if res[-3:] == 'POS' else 0.0 for res in gr_truth.iloc[:, 0] ]
    gr_truth.iloc[:, 0] = test
    gr_truth.columns = ["Truth"]
    return gr_truth

def prepare_data(norm, binary_matrix):
    """ 
    For each column of the expression matrix, multiply the 
    smaller binary matrix by the common genes to get a matrix that is
    basically the same as the binary matrix but with expression instead of
    1's. This leads to a list of matrices.
    """
    list_inputs = []
    bm_tens = torch.tensor(binary_matrix.T.values)
    for sample in norm.columns: 
        item = torch.mul(bm_tens, torch.tensor(norm.loc[common_genes, sample].values))
        list_inputs.append(np.array(item).astype(np.float32))
    return list_inputs 