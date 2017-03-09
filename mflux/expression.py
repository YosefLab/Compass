from __future__ import print_function, division
import numpy as np
import pandas as pd


def load_file(file_name):
    # type: (str) -> SummarizedExperiment
    """
    Loads an expression matrix from file

    File is assumed to be tab-delimited.  Genes on rows, samples on columns.
    Genes and sample labels should be in the file

    Genes normalized to upper-case
    Duplicates are handled by summing values

    If data is not already on a log-range, then it is log2(x+1) transformed
    """

    data = pd.read_table(file_name, index_col=0)
    data.index = data.index.str.upper()


    data_range = data.values.max() - data.values.min()
    if data_range < 30:
        data_is_logged = True
    else:
        data_is_logged = False
    
    # aggregate on genes to eliminate duplicates in case they exist
    # This needs to be done on un-logged data
    if data_is_logged:
        data = 2**data-1
        data_is_logged = False

    data['genekeycolumn'] = data.index
    data = data.groupby('genekeycolumn').sum()

    if not data_is_logged:
        data = np.log2(data + 1)

    se = SummarizedExperiment(data)
    return se


class SummarizedExperiment(object):

    def __init__(self, expression_data, row_data=None, meta_data=None, col_data=None):

        self.expression_data = expression_data
        self.row_data = row_data
        self.col_data = col_data
        self.meta_data = meta_data

        if self.row_data is not None:
            assert self.row_data.shape[0] == self.expression_data.shape[0]

        if self.col_data is not None:
            assert self.col_data.shape[0] == self.expression_data.shape[1]
            