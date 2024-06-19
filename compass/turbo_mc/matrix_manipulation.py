import numpy as np
import pandas as pd

from typing import List, Set, Tuple

NONCONSTANT_COLUMN_THRESHOLD = 1e-3


def are_equal(u: np.array, v: np.array) -> bool:
    r"""
    The business logic for determining whether two columns are equal lives here.
    """
    return np.max(np.abs(u - v)) < NONCONSTANT_COLUMN_THRESHOLD


def is_nonconstant(v: np.array) -> bool:
    r"""
    The business logic for determining whether a feature is nonconstant lives here.
    :param v: 1D array
    :return: True if the column is considered nonconstant.
    """
    return np.abs(np.nanmax(v) - np.nanmin(v)) > NONCONSTANT_COLUMN_THRESHOLD


def is_constant(v: np.array) -> bool:
    r"""
    The business logic for determining whether a feature is constant lives here.
    :param v: 1D array
    :return: True if the column is considered constant.
    """
    return not is_nonconstant(v)


def calc_is_nonconstant_column(X: np.array) -> np.array:
    r"""
    The business logic for determining whether a feature is nonconstant lives here.
    :param X: 2D array
    :return is_nonconstant_column: Boolean np.array with True/False values, one for each column of X.
    """
    return np.abs(np.nanmax(X, axis=0) - np.nanmin(X, axis=0)) > NONCONSTANT_COLUMN_THRESHOLD


def calc_is_constant_column(X: np.array) -> np.array:
    r"""
    The business logic for determining whether a feature is constant lives here.
    :param X: 2D array
    :return is_constant_column: Boolean np.array with True/False values, one for each column of X.
    """
    return ~calc_is_nonconstant_column(X)


def get_nonconstant_columns(X: np.array) -> np.array:
    r"""
    The business logic for determining which features are nonconstant lives here.
    :param X: 2D array
    :return nonconstant_columns: Indices of columns of X that are nonconstant.
    """
    nonconstant_columns = np.where(calc_is_nonconstant_column(X))[0]
    return nonconstant_columns


def get_constant_columns(X: np.array) -> np.array:
    r"""
    The business logic for determining which features are constant lives here.
    :param X: 2D array
    :return constant_columns: Indices of columns of X that are constant.
    """
    constant_columns = np.where(calc_is_constant_column(X))[0]
    return constant_columns


def get_list_of_random_matrix_indices_fast(
        nrows: int,
        ncols: int,
        sampling_density: float = 0.1,
        verbose: bool = False) \
        -> List[Tuple[int, int]]:
    r"""
    Return a random set of indices indexing a matrix of size nrows X ncols.
    Each matrix entry is chosen with probability 'sampling_density'.
    Therefore, some rows and columns might be sampled to different depths.
    If you want to sample each row (resp column) to equal depth, use
    get_list_of_random_matrix_indices instead.
    :param nrows: number of rows of matrix
    :param ncols: number of columns of matrix
    :param sampling_density: number of matrix entries to sample
    :return res: list of indices (r, c)
    """
    rows, cols = np.where(np.random.binomial(1, sampling_density, size=(nrows, ncols)))
    return list(zip(rows, cols))


def get_list_of_random_matrix_indices(
        nrows: int,
        ncols: int,
        sampling_density: float = 0.1,
        verbose: bool = False,
        randomize_rows_and_columns: bool = False) \
        -> List[Tuple[int, int]]:
    r"""
    TODO: cythonize.
    Return a random set of indices indexing a matrix of size nrows X ncols.
    Each row and column is sampled to approximately the same depth 'sampling_density'.
    :param nrows: number of rows of matrix
    :param ncols: number of columns of matrix
    :param sampling_density: number of matrix entries to sample
    :param randomize_rows_and_columns: If True, rows and columns will be shuffled.
        Recommended to set to True (is False by default for backwards compatibility).
    :return res: list of indices (r, c)
    """
    total_entries_sampled = int(nrows * ncols * sampling_density)
    total_entries_sampled_per_column = int(total_entries_sampled / ncols)
    if total_entries_sampled_per_column == 0:
        raise ValueError("No entries would be sampled")
    row_indices = []  # type: List[int]
    while len(row_indices) < 2 * total_entries_sampled:
        row_indices += list(np.random.choice(range(nrows), size=nrows, replace=False))
    key = 0
    res = []  # type: List[Tuple[int, int]]
    hits_per_row = [0] * nrows
    hits_per_col = [0] * ncols
    for c in range(ncols):
        entries_for_this_column = set()  # type: Set[Tuple[int, int]]
        for i in range(total_entries_sampled_per_column):
            while (row_indices[key], c) in entries_for_this_column:
                key += 1
            entries_for_this_column.add((row_indices[key], c))
            hits_per_row[row_indices[key]] += 1
            hits_per_col[c] += 1
            key += 1
        res += list(entries_for_this_column)
    # Check that there are no repetitions
    assert(len(res) == len(set(res)))
    # Check that all columns are equally represented
    assert((max(hits_per_col) - min(hits_per_col)) == 0)
    if verbose:
        print(f"max sampling dicrepancy between rows:\
              {(max(hits_per_row) - min(hits_per_row)) / (min(hits_per_row) + 1e-8)}")
    if randomize_rows_and_columns:
        row_mapping = np.random.choice(range(nrows), nrows, False)
        col_mapping = np.random.choice(range(ncols), ncols, False)
        res = [(row_mapping[r], col_mapping[c]) for (r, c) in res]
    return res


def observe_entries_and_hide_rest(
        X: np.array,
        matrix_indices: List[Tuple[int, int]],
        verbose: bool = False)\
        -> np.array:
    r"""
    Given a matrix X, and a list of indices, returns a matrix with
    matrix_indices observed. Everything else is np.nan
    :param X: pd.DataFrame
    :param matrix_indices: indices of observed entries
    :return res: pd.DataFrame of triples (row, col, val)
    """
    X_observed = np.empty(shape=X.shape)
    X_observed[:] = np.nan
    rows, cols = zip(*matrix_indices)
    X_observed[rows, cols] = X[rows, cols]
    return X_observed


def observe_entries(
        X: np.array,
        matrix_indices: List[Tuple[int, int]],
        verbose: bool = False)\
        -> pd.DataFrame:
    r"""
    Given a matrix X, and a list of indices, returns a DataFrame with
    those (row, col, val) tuples
    :param X: pd.DataFrame
    :param matrix_indices: indices of observed entries
    :return res: pd.DataFrame of triples (row, col, val)
    """
    if verbose:
        print("In observe_entries ... ")  # pragma: no cover
    assert(len(X.shape) == 2)
    rows, cols = zip(*matrix_indices)
    vals = X[rows, cols]
    observations = \
        pd.DataFrame({'row': rows,
                      'col': cols,
                      'val': vals})
    if verbose:
        print("Out observe_entries")  # pragma: no cover
    return observations


def observation_list_to_df(
        observations_list: List[Tuple[int, int, float]])\
        -> pd.DataFrame:
    r"""
    Convert list of (row, col, val) tuples to a DataFrame.
    Observations are sorted by column first, then row.
    :param observations: list of (row, col, val)
    :return res: pd.DataFrame of triples (row, col, val)
    """
    # Check that observations are all distinct
    assert(len(set([(r, c) for (r, c, _) in observations_list]))
           == len(observations_list))
    observations = \
        pd.DataFrame({'row': [r for r, _, _ in observations_list],
                      'col': [c for _, c, _ in observations_list],
                      'val': [v for _, _, v in observations_list]})
    observations.sort_values(by=['col', 'row'], inplace=True, ignore_index=True)
    return observations


def matrix_from_observations(
        observations: pd.DataFrame,
        nrows: int,
        ncols: int) -> np.array:
    r"""
    Given a set of observations, return the matrix. pd.nan is used for unobserved entries.
    The number of rows and columns of the resulting matrix is required. Although
    they could be guessed from 'observations', I want to be pedantic here to avoid
    bugs where the last row of the matrix was completely unobserved, which might be
    the case e.g. with a linear regression model learnt from only a few rows.
    :param observations: pd.DataFrame of (rol, col, val) triples
    :param nrows: Number of rows of the matrix.
    :param ncols: Number of columns of the matrix.
    :return res: np.array with np.nan where the matrix is not observed.
    """
    res = np.zeros(shape=(nrows, ncols))
    res[:, :] = np.nan
    res[observations['row'], observations['col']] = observations['val']
    return res


def get_observations_from_partially_observed_matrix(X: np.array)\
        -> List[Tuple[int, int, float]]:
    observed_indices = np.where(~np.isnan(X))
    rows = list(observed_indices[0])
    cols = list(observed_indices[1])
    vals = list(X[observed_indices])
    return observation_list_to_df(list(zip(rows, cols, vals)))


def complete_matrix(
        X_observed: np.array,
        model,  # MatrixCompletionModel; not typing it to avoid circular import...
        verbose: bool = False)\
        -> None:
    r"""
    TODO: This is super slow... slower than training the model!
    Given a partially observed matrix X_observed and a MatrixCompletionModel,
    returns a new matrix where pd.nans in X_observed are imputed with the model.
    :param X_observed: Partially observed matrix
    :param model: A MatrixCompletionModel
    :return X_completion: Completing of X_observed based on model
    """
    if verbose:
        print("In complete_matrix ... ")  # pragma: no cover
    nrows, ncols = X_observed.shape
    X_completion = X_observed.copy()
    for r in range(nrows):
        for c in range(ncols):
            if pd.isna(X_completion[r, c]):
                X_completion[r, c] = model.predict(r, c)
    if verbose:
        print("Out complete_matrix")  # pragma: no cover
    return X_completion


def standardize(X: np.array) -> np.array:
    r"""
    Constant columns are standardized to 0. Remaining columns are standardized by
    subtracting the mean and dividing by the standard deviation.
    :param X: 2D np.array of size observations X features.
    :param X_normalized: (X - column_means) / column_stds
    """
    X = X.copy()
    nonconstant_columns = get_nonconstant_columns(X)
    X_nonconstant = X[:, nonconstant_columns]
    column_means = np.expand_dims(np.mean(X_nonconstant, axis=0).flatten(), axis=0)
    column_stds = np.expand_dims(np.std(X_nonconstant, axis=0).flatten(), axis=0)
    X_nonconstant_standardized = (X_nonconstant - column_means) / column_stds
    X_standardized = np.zeros_like(X)
    X_standardized[:, nonconstant_columns] = X_nonconstant_standardized
    return X_standardized
