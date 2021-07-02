"""
utils.py

Functions to be used by other optimization routines
"""

from __future__ import print_function, division
import cplex
import scipy.io
import pandas as pd

def get_steadystate_constraints(model):
    """
    Uses the s_mat to define connectivity constraints
    """
    s_mat = model.getSMAT()

    lin_expr = []
    senses = []
    rhs = []
    names = []

    for metab, rx in s_mat.items():
        if len(rx) == 0:
            continue

        ind = [x[0] for x in rx]
        val = [x[1] for x in rx]

        lin_expr.append(cplex.SparsePair(ind=ind, val=val))
        senses.append('E')
        rhs.append(0)
        names.append(metab)

    return lin_expr, senses, rhs, names


def reset_objective(problem):
    """
    Clears all the objective coefficients for the current problem
    by setting all to 0
    """

    names = problem.variables.get_names()
    zeros = [0 for x in names]

    problem.objective.set_linear(zip(names, zeros))
    problem.objective.set_name('none')

def read_data(data):
    if len(data) == 1:
        return pd.read_csv(data[0], sep='\t', index_col=0)
    else:
        return read_mtx(data[0], data[1], data[2])

def read_sample_names(data):
    if len(data) == 1:
        return pd.read_csv(data[0], sep='\t', index_col=0, nrows=1).columns
    elif len(data) >= 3 and data[2] is not None:
        return pd.read_csv(data[2], sep='\t', header=None).to_numpy().ravel()
    else:
        #Sample names not provided
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

def read_knn(knn_data):
    """
        Reads a knn_file in either csv or sparse format
    """
    if len(knn_data) > 2:
        return read_mtx(knn_data[0], knn_data[1], knn_data[2])
    elif len(knn_data) > 1:
        return read_mtx(knn_data[0], knn_data[1], knn_data[1])
    else:
        return pd.read_csv(knn_data[0], sep='\t', index_col=0)