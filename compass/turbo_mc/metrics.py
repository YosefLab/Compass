r"""
The basic API of a metric is (X, X_observed, X_pred, exclude_observed)
"""
import numpy as np
import scipy.stats
# import scipy.stats  # I could compute r2 with scipy.stats.linregress but rather do it by hand.
from scipy.sparse.linalg import svds
from typing import Optional

from .matrix_manipulation import standardize, are_equal, is_constant


def compute_r2(
        y_true: np.array,
        y_observed: np.array,
        y_pred: np.array,
        exclude_observed: bool = True) -> float:
    r"""
    The R2 is computed a 1.0 - ss_res / (ss_tot + 1e-16) where ss_res is the residual sum of
    squares, and ss_tot is the sum of squares of the baseline model np.nanmean(y_observed).
    The R2 is clipped -1.0
    Note that, as for MSE computation:
        - predicted columns that are almost equal to the ground truth are declared R2 1.0
        - predicted columns that are not almost equal yet ground truth is constant, have R2 -1.0.
        - for all other columns, the usual R2 is computed
    """
    baseline_mean = np.nanmean(y_observed)
    # _, _, r2, _, _ = scipy.stats.linregress(y_true, y_pred)
    indices_to_use = np.where(np.isnan(y_observed)) if exclude_observed else\
        np.where(~np.isnan(y_true))  # Hacky way to get all indices...
    # Border case 0: all the column entries are observed so there are no errors.
    if len(indices_to_use[0]) == 0:
        return 1.0
    # Border case 1: predictions are very close to ground truth
    if are_equal(y_true[indices_to_use], y_pred[indices_to_use]):
        return 1.0
    # Border case 2: predictions are not very close to ground truth, yet ground truth is constant
    if not are_equal(y_true[indices_to_use], y_pred[indices_to_use]) and is_constant(y_true):
        return -1.0
    # Usual computation:
    ss_tot = np.sum((y_true[indices_to_use] - baseline_mean) ** 2)
    ss_res = np.sum((y_true[indices_to_use] - y_pred[indices_to_use]) ** 2)
    r2 = 1.0 - ss_res / ss_tot
    r2 = max(r2, -1.0)
    return r2


def compute_r2s(
        X: np.array,
        X_observed: np.array,
        X_completion: np.array,
        exclude_observed: bool = True,
        verbose: bool = False)\
        -> np.array:
    r"""
    Computes r2 accross all columns, only on UN-observed entries. Truncated below by -1.0.
    The only reason I care about considering observed entries for MSE is the Theory behind it.
    Here there is no theory so I won't even give the option of computing R2 with observed entries.
    """
    if verbose:
        print(f"In compute_r2s ... ")  # pragma: no cover
    assert(X.shape == X_observed.shape)
    assert(X.shape == X_completion.shape)
    nrows, ncols = X.shape
    r2s = np.zeros(shape=(ncols))
    for c in range(ncols):
        y_true = X[:, c]
        y_observed = X_observed[:, c]
        y_pred = X_completion[:, c]
        r2 = compute_r2(y_true, y_observed, y_pred, exclude_observed=exclude_observed)
        r2s[c] = r2
    if verbose:
        print(f"Out compute_r2s")  # pragma: no cover
    return r2s


def compute_pearson_r2(
        y_true: np.array,
        y_observed: np.array,
        y_pred: np.array,
        exclude_observed: bool = True) -> float:
    r"""
    Standard Pearson R2 correlation: a number between -1 and 1.
    Note that:
        - predicted columns that are almost equal to the ground truth are declared Pearson R2 1.0
        - predicted columns that are not almost equal yet ground truth is constant, have Pearson R2 0.0.
        - for all other columns, the usual Pearson R2 is computed
    """
    indices_to_use = np.where(np.isnan(y_observed)) if exclude_observed else\
        np.where(~np.isnan(y_true))  # Hacky way to get all indices...
    # Border case 0: all the column entries are observed so there are no errors.
    if len(indices_to_use[0]) == 0:
        return 1.0
    # Border case 1: predictions are very close to ground truth
    if are_equal(y_true[indices_to_use], y_pred[indices_to_use]):
        return 1.0
    # Border case 2: predictions are not very close to ground truth, yet ground truth is constant
    if not are_equal(y_true[indices_to_use], y_pred[indices_to_use]) and is_constant(y_true):
        return 0.0
    # Usual computation:
    # Border case: all true values or predicted values are constant (but the other is not)
    if np.std(y_true[indices_to_use]) == 0.0:
        # assert(np.std(y_pred[indices_to_use]) != 0.0)  # This has blown up in my face, np.std doesn't seem robust.
        return 0.0
    if np.std(y_pred[indices_to_use]) == 0.0:
        # assert(np.std(y_true[indices_to_use]) != 0.0)  # Hasn't blown up in my face, but commenting out.
        return 0.0
    pearson_r2 = np.corrcoef(y_true[indices_to_use], y_pred[indices_to_use])[0, 1]
    # Paranoid check so this doesn't ruin my day. TODO: See if this ever hits.
    if np.isnan(pearson_r2):
        # print("WARNING: np.corrcoef returned np.nan")
        return 0.0
    return pearson_r2


def compute_pearson_r2s(
        X: np.array,
        X_observed: np.array,
        X_completion: np.array,
        exclude_observed: bool = True)\
        -> np.array:
    r"""
    Computes Pearson R2 accross all columns.
    Note that:
        - predicted columns that are almost equal to the ground truth are declared Pearson R2 1.0
        - predicted columns that are not almost equal yet ground truth is constant, have Pearson R2 0.0.
        - for all other columns, the usual Pearson R2 is computed
    """
    assert(X.shape == X_observed.shape)
    assert(X.shape == X_completion.shape)
    nrows, ncols = X.shape
    pearson_r2s = np.zeros(shape=(ncols))
    for c in range(ncols):
        y_true = X[:, c]
        y_observed = X_observed[:, c]
        y_pred = X_completion[:, c]
        pearson_r2 = compute_pearson_r2(y_true, y_observed, y_pred, exclude_observed=exclude_observed)
        pearson_r2s[c] = pearson_r2
    return pearson_r2s


def compute_spearman_r2(
        y_true: np.array,
        y_observed: np.array,
        y_pred: np.array,
        exclude_observed: bool = True) -> float:
    r"""
    Standard Spearman R2 correlation: a number between -1 and 1.
    Note that:
        - predicted columns that are almost equal to the ground truth are declared Spearman R2 1.0
        - predicted columns that are not almost equal yet ground truth is constant, have Spearman R2 0.0.
        - for all other columns, the usual Spearman R2 is computed
    """
    indices_to_use = np.where(np.isnan(y_observed)) if exclude_observed else\
        np.where(~np.isnan(y_true))  # Hacky way to get all indices...
    y_true_sub = y_true[indices_to_use]
    y_pred_sub = y_pred[indices_to_use]
    # Border case 0: all the column entries are observed so there are no errors.
    if len(indices_to_use[0]) == 0:
        return 1.0
    # Border case 1: predictions are very close to ground truth
    if are_equal(y_true_sub, y_pred_sub):
        return 1.0
    # Border case 2: predictions are not very close to ground truth, yet ground truth is constant
    if not are_equal(y_true_sub, y_pred_sub) and is_constant(y_true):
        return 0.0
    # Usual computation:
    # Border case: all true values or predicted values are constant (but the other is not)
    if np.std(y_true_sub) == 0.0:
        # assert(np.std(y_pred_sub) != 0.0)  # This has blown up in my face, np.std doesn't seem very robust.
        return 0.0
    if np.std(y_pred_sub) == 0.0:
        # assert(np.std(y_true_sub) != 0.0)  # Hasn't blown up in my face, but commenting out.
        return 0.0
    spearman_r2 = scipy.stats.spearmanr(y_true_sub, y_pred_sub).correlation
    # Even then, sometimes spearman_r2 can return np.nan :( Not sure why
    if np.isnan(spearman_r2):
        # print("WARNING: scipy.stats.spearmanr returned np.nan")
        return 0.0
    return spearman_r2


def compute_spearman_r2s(
        X: np.array,
        X_observed: np.array,
        X_completion: np.array,
        exclude_observed: bool = True)\
        -> np.array:
    r"""
    Computes Spearman R2 accross all columns.
    Note that:
        - predicted columns that are almost equal to the ground truth are declared Spearman R2 1.0
        - predicted columns that are not almost equal yet ground truth is constant, have Spearman R2 0.0.
        - for all other columns, the usual Spearman R2 is computed
    """
    assert(X.shape == X_observed.shape)
    assert(X.shape == X_completion.shape)
    nrows, ncols = X.shape
    spearman_r2s = np.zeros(shape=(ncols))
    for c in range(ncols):
        y_true = X[:, c]
        y_observed = X_observed[:, c]
        y_pred = X_completion[:, c]
        spearman_r2 = compute_spearman_r2(y_true, y_observed, y_pred, exclude_observed=exclude_observed)
        spearman_r2s[c] = spearman_r2
    return spearman_r2s


def compute_mse(
        y_true: np.array,
        y_observed: np.array,
        y_pred: np.array,
        give_each_column_equal_weight: bool = False,
        exclude_observed: bool = True,
        return_sum_instead_of_mean: bool = False) -> float:
    r"""
    Given true values, observed values, and predicted values, compute the MSE.
    :param give_each_column_equal_weight: If True:
        - predicted columns that are almost equal to the ground truth are declared MSE 0.
        - predicted columns that are not almost equal yet ground truth is constant, have MSE nan.
        - for all other columns, the usual MSE is scaled by dividing by the true std ** 2.
    :param exclude_observed: If True, observed values will be excluded from the calculation.
        Inclusion is usefull for comparing against the Theoretical minimum low-rank
        reconstruction error with SVD.
    :param return_sum_instead_of_mean: If True, will return sum of errors instead of the mean.
    """
    res = 0.0
    indices_to_use = np.where(np.isnan(y_observed)) if exclude_observed\
        else np.where(~np.isnan(y_true))  # hacky way to say "all indices"...
    # Border case 0: all the column are observed so there are no errors.
    if len(indices_to_use[0]) == 0:
        return 0.0
    error = y_true[indices_to_use] - y_pred[indices_to_use]
    res = np.mean(error * error)
    if give_each_column_equal_weight:
        # Border case 1: predictions are very close to ground truth
        if are_equal(y_true[indices_to_use], y_pred[indices_to_use]):
            return 0.0
        # Border case 2: predictions are not very close to ground truth, yet ground truth is constant
        if not are_equal(y_true[indices_to_use], y_pred[indices_to_use]) and is_constant(y_true):
            return np.nan
        res /= np.std(y_true) ** 2
    if return_sum_instead_of_mean:
        res *= len(error)
    return res


def compute_mse_per_column(
        X: np.array,
        X_observed: np.array,
        X_completion: np.array,
        exclude_observed: bool = True,
        give_each_column_equal_weight: bool = False,
        verbose: bool = False,
        return_sum_instead_of_mean: bool = False)\
        -> np.array:
    r"""
    Main method for computing MSEs in various forms (e.g. variance explained,
    variance explained excluding observed)
    For each column, computes the MSE between predictions and ground truth.
    :param exclude_observed: If True, observed values will be excluded from the calculation.
        Inclusion is usefull for comparing against the Theoretical minimum low-rank
        reconstruction error with SVD.
    :param give_each_column_equal_weight: If True:
        - predicted columns that are almost equal to the ground truth are declared MSE 0.
        - predicted columns that are not almost equal yet ground truth is constant, have MSE nan.
        - for all other columns, the usual MSE is scaled by dividing by the true std ** 2.
    :param return_sum_instead_of_mean: If True, will return sum of errors instead of the mean.
    """
    ncols = X.shape[1]
    res = np.zeros(shape=(ncols))
    for c in range(ncols):
        y_true = X[:, c]
        y_observed = X_observed[:, c]
        y_pred = X_completion[:, c]
        res[c] = compute_mse(
            y_true,
            y_observed,
            y_pred,
            give_each_column_equal_weight,
            exclude_observed,
            return_sum_instead_of_mean)
    return res


def compute_mse_giving_equal_weight_to_all_columns(
        X: np.array,
        X_observed: np.array,
        X_completion: np.array,
        exclude_observed: bool = True,
        verbose: bool = False)\
        -> float:
    r"""
    TODO: Remove?
    Computes the MSE, giving the same weight to each column. This is achieved by dividing each
    column by the true standard deviation.
    """
    mse_per_column = compute_mse_per_column(
        X=X,
        X_observed=X_observed,
        X_completion=X_completion,
        exclude_observed=exclude_observed,
        give_each_column_equal_weight=True,
        verbose=verbose
    )
    ncols = X.shape[1]
    mse_giving_equal_weight_to_all_columns = 0.0
    total_entries_summed = 0
    for c in range(ncols):
        indices_to_use = np.where(np.isnan(X_observed[:, c])) if exclude_observed\
            else np.where(~np.isnan(X[:, c]))  # hacky way to say "all indices"...
        mse_giving_equal_weight_to_all_columns += mse_per_column[c] * len(indices_to_use[0])
        total_entries_summed += len(indices_to_use[0])
    if total_entries_summed > 0:
        mse_giving_equal_weight_to_all_columns /= total_entries_summed
    return mse_giving_equal_weight_to_all_columns


def compute_mse_global(
        X: np.array,
        X_observed: np.array,
        X_completion: np.array,
        exclude_observed: bool = True,
        verbose: bool = False)\
        -> float:
    r"""
    Computes the global MSE between predictions and ground truth.
    Note that this is not per-column: it is the _average_ over all columns.
    """
    if verbose:
        print(f"In compute_mse_global ... ")  # pragma: no cover
    indices_to_use = np.where(np.isnan(X_observed)) if exclude_observed\
        else np.where(~np.isnan(X))  # hacky way to say "all indices"...
    y_true = X[indices_to_use]
    y_pred = X_completion[indices_to_use]
    error = y_true - y_pred
    if verbose:
        print(f"Out compute_mse_global")  # pragma: no cover
    if len(error) == 0:  # Border case
        return 0.0
    return np.mean(error * error)


def variance_explained_upper_bounds(
        X: np.array,
        give_each_column_equal_weight: bool = True,
        k: Optional[int] = None) -> np.array:
    r"""
    TODO: At position 0, make the variance explained be zero?
    Maximum possible variance explained when approximating X by a rank k matrix,
    for all k. Constant columns are replaced by zero if give_each_column_equal_weight=True.
    :param give_each_column_equal_weight: If True, constant columns are dropped,
        and columns are divided by their standard deviation.
    :param k: If specified, the truncated SVD will be computed. Useful for large datasets where
        a full SVD would not be possible.
    :return: 1D np.array which, at position k-1, contains the maximum possible
        variance explained by a rank k approximation.
    """
    if give_each_column_equal_weight:
        X = standardize(X)
    if k is not None:
        U, Sigma, V = svds(X, k=k)
        Sigma = Sigma[::-1]
    else:
        U, Sigma, V = np.linalg.svd(X, full_matrices=False)
    variances_explained = np.cumsum(Sigma * Sigma) / np.sum(Sigma * Sigma)
    return variances_explained


def compute_best_low_rank_approximation(X: np.array, k: int) -> np.array:
    r"""
    Computes the best rank k approximation to X (as given by either Frobenius or
    operator norm). It consists of just keeping the top k singular values of X.
    :param X: 2D np.array matrix to approximate.
    :param k: Rank of the low-rank approximation.
    :return: 2D np.array, the best rank k approximation to X.
    """
    U, Sigma, V = np.linalg.svd(X, full_matrices=False)
    res = U[:, :k] @ np.diag(Sigma[:k]) @ V[:k, :]
    return res


def compute_variance_explained(
        X: np.array,
        X_pred: np.array,
        give_each_column_equal_weight: bool = True) -> float:
    r"""
    Implemented as sum of MSEs of all columns, divided by squared L2 norm of X.
    :param give_each_column_equal_weight: If True:
        - predicted columns that are almost equal to the ground truth are declared MSE 0.
        - predicted columns that are not almost equal yet ground truth is constant, have MSE nan.
        - for all other columns, the usual MSE is scaled by dividing by the true std ** 2.
    """
    sses = compute_mse_per_column(X, X, X_pred, exclude_observed=False,
                                  give_each_column_equal_weight=give_each_column_equal_weight,
                                  return_sum_instead_of_mean=True)
    if give_each_column_equal_weight:
        X = standardize(X)
    res = 1.0 - np.sum(sses) / np.sum(X ** 2)
    return res


def compute_variance_explained_excluding_observed_entries(
        X: np.array,
        X_observed: np.array,
        X_pred: np.array,
        give_each_column_equal_weight: bool = True) -> float:
    r"""
    TODO: Untested & Unused.
    Implemented as sum of MSEs of all columns, divided by squared L2 norm of X.
    :param give_each_column_equal_weight: If True:
        - predicted columns that are almost equal to the ground truth are declared MSE 0.
        - predicted columns that are not almost equal yet ground truth is constant, have MSE nan.
        - for all other columns, the usual MSE is scaled by dividing by the true std ** 2.
    """
    sses = compute_mse_per_column(X, X_observed, X_pred, exclude_observed=True,
                                  give_each_column_equal_weight=give_each_column_equal_weight,
                                  return_sum_instead_of_mean=True)
    if give_each_column_equal_weight:
        X = standardize(X)
    indices_to_use = np.where(np.isnan(X_observed))
    if len(indices_to_use[0]) == 0:  # Border case: nothing to do
        return 1.0
    res = 1.0 - np.sum(sses) / np.sum(X[indices_to_use] ** 2)
    return res


def compute_squared_error_fast(
        X_true: np.array,
        X_completion: np.array) -> float:
    r"""
    Squared Error, normalizing all columns. X_true and X_completion must be fully observed.
    Can blow up if constant functions are not exactly constant in the predictions. In that case, use the
    slower compute_mse_giving_equal_weight_to_all_columns.
    """
    col_means = X_true.mean(axis=0)
    col_stds = X_true.std(axis=0)
    X_true = (X_true - col_means) / (col_stds + 1e-16)
    X_completion = (X_completion - col_means) / (col_stds + 1e-16)
    return np.linalg.norm(X_true - X_completion, 'fro') ** 2


def compute_variance_explained_fast(
        X_true: np.array,
        X_completion: np.array,
        pseudo_std: float = 1e-3) -> float:
    r"""
    Variance Explained, normalizing all columns. X_true and X_completion must be fully observed.
    Can blow up if constant functions are not exactly constant in the predictions. In that case, use the
    slower compute_variance_explained.
    :param pseudo_std: Added to the std when normalizing columns by their std.
        Used to avoid producing np.nans (infinities) in border cases where
        numerical instability leads to standard deviations of 0 that are really
        infinitesimally small but not zero.
    """
    col_means = X_true.mean(axis=0)
    col_stds = X_true.std(axis=0)
    X_true = (X_true - col_means) / (col_stds + pseudo_std)
    X_completion = (X_completion - col_means) / (col_stds + pseudo_std)
    return 1.0 - np.linalg.norm(X_true - X_completion, 'fro') ** 2 / (np.linalg.norm(X_true, 'fro') ** 2 + pseudo_std)
