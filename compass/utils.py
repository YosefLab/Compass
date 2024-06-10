"""
utils.py

Functions to be used by other optimization routines
"""

from __future__ import print_function, division
import scipy.io
import pandas as pd
import numpy as np
from .globals import MODEL_DIR
import os
import anndata

import gurobipy as gp

def get_steadystate_constraints(model, gp_model):

    """
    Uses the s_mat to define connectivity constraints
    """

    # s_mat is a dictionary with metabolites as keys and list of (reactions, stoichiometric coefficients) as values
    s_mat = model.getSMAT()

    lin_expr = []
    rhs = []
    names = []

    for metab, rx in s_mat.items():
        # If there is no reaction associated with the given metabolite, then skip
        if len(rx) == 0:
            continue

        # x[0] is name of reaction
        # x[1] is stoichiometric coefficient of metabolite in reaction x[0]
        expr = gp.LinExpr()
        for x in rx:
            expr += x[1] * gp_model.getVarByName(x[0])

        lin_expr.append(expr)
        rhs.append(0)
        names.append(metab)

    return lin_expr, rhs, names


def reset_objective(gp_model):
    """
    Clears all the objective coefficients for the current problem
    by setting all to 0
    """

    gp_model.setObjective(0)
    gp_model.update()
    
def read_data(data):
    if len(data) == 1:
        ext = os.path.splitext(data[0])[-1]
        if ext == '.h5ad':
            return anndata.read_h5ad(data[0]).to_df().T
        else:
            return pd.read_csv(data[0], sep='\t', index_col=0)
    else:
        return read_mtx(data[0], data[1], data[2])

def read_annotations(data):
    if len(data) == 1:
        ext = os.path.splitext(data[0])[-1]
        if ext == '.h5ad':
            res = anndata.read_h5ad(data[0])[:,:0].copy() #Slice to remove all gene observations
            return res
        else:
            return None
    elif len(data) == 3:
        return None

def write_output(output, path, args):
    if args['anndata_output']:
        #TODO: Add more control over output format
        
        #Output will be indexed by "sample_%d".format(index) unless reading in slow names
        #output.var = args['anndata_annotations'].obs

        #Generally only observational annotations are relevant after Compass algorithm
        annot = args['anndata_annotations']
        res = anndata.AnnData(X=output.T, obs=annot.obs, uns=annot.uns, obsm=annot.obsm, obsp=annot.obsp)
        res.write(path+'.h5ad', compression='gzip') 
    else:
        output.to_csv(path+".tsv", sep="\t")

def read_sample_names(data, slow_names=True):
    """
    Reads in sample names for dataset

    Some data input formats do not support fast ways to read sample names (h5ad) and when slow_names is False, reading them will be skipped.
    """
    if len(data) == 1:
        ext = os.path.splitext(data[0])[-1]
        if ext == '.h5ad':
            if slow_names:
                return anndata.read_h5ad(data[0]).obs.index
        else:
            return pd.read_csv(data[0], sep='\t', index_col=0, nrows=1).columns
    elif len(data) >= 3 and data[2] is not None:
        return pd.read_csv(data[2], sep='\t', header=None).to_numpy().ravel()
    #Sample names not provided or not efficient to read
    return None

def indexed_sample_names(n):
    return ['sample_'+str(i) for i in range(n)]

def read_mtx(mtx_file, rows_file, columns_file=None):
    """
        Reads an mtx file into a pandas dataframe for doing Compass stuff with. Primarily for reading gene expression files.
    """
    mtx = scipy.io.mmread(mtx_file)
    rows = pd.read_csv(rows_file, sep='\t', header=None)
    if columns_file is not None:
        columns = pd.read_csv(columns_file, sep='\t', header=None).to_numpy().ravel()
    else:
        columns = indexed_sample_names(mtx.shape[1])
    if pd.__version__ >= '1':
        return pd.DataFrame.sparse.from_spmatrix(mtx, index=rows.to_numpy().ravel(), columns = columns)
    else:
        return pd.SparseDataFrame(mtx, index=rows.to_numpy().ravel(), columns = columns)

def read_knn(knn_data, data=None, dist=False):
    """
        Parses a knn input format from sklearn's nearest neighbors format (possibly wrapped in a Pandas dataframe)
        Returns None if the result does not match
    """
    if (knn_data.endswith("npy")):
        knn = np.load(knn_data)
        if knn.shape[0] == data.shape[1]:
            return knn
        else:
            return None
    else:
        # knn should be of shape (# of samples, k+1)
        # first column is the sample names of the k nearest neighbors
        # following k columns are indices of the k nearest neighbors
        # data is of shape (# of genes, # of cells)
        knn = pd.read_csv(knn_data, sep='\t', index_col=0)
        if data is None:
            return knn.values #No choice but to assume that the indices are the same as input data
        elif len(knn.index) != len(data.columns):
            return None
        elif np.all(knn.index == data.columns):
            return knn.values 
        elif np.all(np.sort(knn.index) == np.sort(data.columns)):
            if dist:
                return knn.loc[data.columns].values
            else:
                #Need this only for the array of indices (because they may be permuted)
                LUT = {}
                for i in range(data.shape[1]):
                    LUT[i] = data.columns.get_loc(knn.index[i])
                return knn.applymap(lambda x: LUT[x]).loc[data.columns].values 
        else:
            return None

def read_knn_ind(knn_data, data=None):
    return read_knn(knn_data, data, dist=False)

def read_knn_dist(knn_data, data=None):
    return read_knn(knn_data, data, dist=True)

def read_metadata(model_name):
    top_dir = os.path.join(MODEL_DIR, model_name)
    metadata_dir = os.path.join(top_dir, 'metadata')
    reaction_metadata_path = os.path.join(metadata_dir, 'reaction_metadata.csv')

    return pd.read_csv(reaction_metadata_path, index_col=0)

def parse_gurobi_license_file(file_path):
    credentials = {}

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("WLSACCESSID="):
                credentials['WLSACCESSID'] = line.split('=')[1]
            elif line.startswith("WLSSECRET="):
                credentials['WLSSECRET'] = line.split('=')[1]
            elif line.startswith("LICENSEID="):
                credentials['LICENSEID'] = int(line.split('=')[1])

    return credentials
