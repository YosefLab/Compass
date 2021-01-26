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

def read_mtx(mtx_file, genes_file, barcodes_file=None):
    """
        Reads an mtx file into a pandas dataframe for doing Compass stuff with
    """
    mtx = scipy.io.mmread(mtx_file)
    genes = pd.read_csv(genes_file, sep='\t', header=None)
    if barcodes_file is not None:
        barcodes = pd.read_csv(barcodes_file, sep='\t', header=None)
        if pd.__version__ >= '1':
            return pd.DataFrame.sparse.from_spmatrix(mtx, index=genes.to_numpy().ravel(), columns = barcodes.to_numpy().ravel())
        else:
            return pd.SparseDataFrame(mtx, index=genes.to_numpy().ravel(), columns = barcodes.to_numpy().ravel())
    else:
        if pd.__version__ >= '1':
            return pd.DataFrame.sparse.from_spmatrix(mtx, index=genes.to_numpy().ravel())
        else:
            return pd.SparseDataFrame(mtx, index=genes.to_numpy().ravel())

        