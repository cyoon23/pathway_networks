import pandas as pd
import numpy as np
import h5py
import urllib.request


def prediction_pipeline(h5_file, libraries, cor):
    """
    Loops through desired libraries from Enrichr, gets pathway-gene libraries,
    creates binary matrices, computes correlation and predicts gene-pathway 
    correlation for new genes of inputted correlation matrix.
    """
    h = h5py.File(h5_file, "a")
    for lib in libraries:
        function_to_genes, gene_set = gene_set_dictionaries(lib)
        bm = gs_binary_matrix(function_to_genes, gene_set)
        gslib = gene_set_library(bm, gene_set)
        pm = prediction_matrix(cor, gslib)
        h.create_dataset("cor_"+lib, data=pm)
        index = (pd.DataFrame(cor.index)).values.astype("S").tolist()
        col = (pd.DataFrame(gslib.columns)).values.astype("S").tolist()
        h.create_dataset("genes_"+lib, data=index)
        h.create_dataset("pathways_"+lib, data=col)
    h.close()


def gene_set_dictionaries(library):
    """
    Extracts gene set library from Enrichr and makes dictionary with functions/pathways
    as keys and genes as values.
    """
    enrichr_url = 'https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName='
    data = urllib.request.urlopen(enrichr_url + library)
    function_to_genes = {}
    gene_set = set() 
    for line in data:
        lst = (str(line.strip())[2:-1]).split(r'\t')
        function = lst[0]
        genes = lst[2:]
        function_to_genes[function] = genes
        gene_set.update(set(genes))
    return function_to_genes, sorted(gene_set)


def gs_binary_matrix(function_to_genes, gene_set):
    """
    Compute binary matrix based on gene's known functions/pathways.
    """
    binary_matrix = np.zeros((len(gene_set), len(function_to_genes)))
    binary_matrix = pd.DataFrame(data=binary_matrix, index=gene_set, columns=list(function_to_genes.keys()))
    for function in binary_matrix.columns: 
        gene_list = function_to_genes[function]
        binary_matrix.loc[gene_list, function] += 1
    return binary_matrix


def gene_set_library(binary_matrix, gene_set):
    """
    Computes correlation matrix for genes of gene set library.
    """
    gslib = np.corrcoef(binary_matrix)
    return pd.DataFrame(gslib, index=gene_set, columns=gene_set)


def prediction_matrix(cor_matrix, gslib):
    """ 
    Computes prediction of gene-pathway correlation.
    """
    idx = [ g for g in gslib.index if g in cor_matrix.index ]
    smaller_cor = cor_matrix.loc[idx] 
    smaller_gslib = gslib.loc[idx]
    pred_matrix = np.dot(np.transpose(smaller_cor), smaller_gslib)
    pred_matrix = pd.DataFrame(pred_matrix, index=smaller_cor.columns, columns=smaller_gslib.columns)
    for gene in smaller_gslib.index:
        row = gslib.loc[gene]
        val = gslib.loc[gene][gene]
        function_sums = np.array(np.sum(row)) - val
        pred_matrix.loc[gene] /= function_sums
    return pd.DataFrame(data=np.transpose(pred_matrix), index=gslib.columns, columns=smaller_gslib.columns)