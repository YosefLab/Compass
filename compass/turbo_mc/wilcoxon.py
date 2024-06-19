import numpy as np
import pandas as pd
import scipy.stats


def wilcoxon_test(x: np.array, y: np.array, return_null_distribution=False):
    '''
    Implements Wilcoxon rank sum test vectorized (as in R).
    Note: The 0.5 'continuity correction':
        https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/wilcox.test.R#L320
        is NOT performed in this method. It is performed in wilcoxon_score
    :param x: vector of values
    :param y: vector of 0/1 respective labels identifying the two groups.
    :param return_null_distribution: if False, only returns correctly_ordered_pairs (useful for fast testing
        this against wilcoxon_test_naive). If True, also returns mu, sigma, and sigma_corrected: the
        parameters for the Normal distribution that approximates the null distribution.
        (return_null_distribution is False by default due to backwards compatibility with my test.)
    :return correctly_ordered_pairs: Number of correctly ordered pairs
    :return mu: mean under null
    :return sigma: std under null, without correction
    :return sigma_corr: std under null
    R implementation: https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/wilcox.test.R#L309
    '''
    assert(len(x) == len(y))
    assert(all((y == 0) | (y == 1)))
    n = len(x)
    n1 = y.sum()
    n0 = n - n1
    correctly_ordered_pairs = (y * scipy.stats.rankdata(x)).sum() - n1 * (n1 + 1) / 2.0
    uniques = np.unique(x, return_counts=True)[1]
    correction = (uniques * (uniques * uniques - 1)).sum()
    assert(n0 + n1 == n)
    mu = n0 * n1 / 2
    sigma = np.sqrt(n0 * n1 * (n + 1) / 12.0)
    sigma_corr = np.sqrt(n0 * n1 / 12.0 * (n + 1 - correction / (n * (n - 1.0))))
    if return_null_distribution:
        return correctly_ordered_pairs, mu, sigma, sigma_corr
    else:
        return correctly_ordered_pairs


def wilcoxon_test_naive(x, y):
    correctly_ordered_pairs = 0
    assert(len(x) == len(y))
    n = len(x)
    for i in range(n):
        for j in range(i + 1, n, 1):
            if x[i] == x[j] and y[i] != y[j]:
                correctly_ordered_pairs += 0.5
            if x[i] < x[j] and y[i] < y[j]:
                correctly_ordered_pairs += 1.0
            if x[i] > x[j] and y[i] > y[j]:
                correctly_ordered_pairs += 1.0
    # n1 = sum(y)
    # n0 = n - n1
    # mu = n0 * n1 / 2
    return correctly_ordered_pairs


def wilcoxon_score(x, y, correct_z_score=True):
    '''
    Returns the wilcoxon rank sum test z score
    :param x: observed values
    :param y: labels
    :param correct_z_score: If True, continuity correction is applied, i.e. 0.5 is added/subtracted from
        the (unnormalized) z score. (prior to division by sigma_corrected) See:
        https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/wilcox.test.R#L320
    '''
    correctly_ordered_pairs, mu, sigma, sigma_corrected =\
        wilcoxon_test(x, y, return_null_distribution=True)
    z_score = (correctly_ordered_pairs - mu) / (sigma_corrected + 1e-16)
    if correct_z_score:
        z_score = z_score - np.sign(z_score) * 0.5 / (sigma_corrected + 1e-16)  # Note that can't overshoot 0.
    return z_score


def calc_cohens_d(x1, x2) -> float:
    r"""
    Cohen's d as per
    https://en.wikipedia.org/wiki/Effect_size#Cohen's_d
    """
    x1, x2 = np.array(x1).flatten(), np.array(x2).flatten()
    n_A, n_B, mu_A, mu_B, var_A, var_B = len(x1), len(x2), x1.mean(), x2.mean(), x1.std() ** 2, x2.std() ** 2
    pooled_sd = np.sqrt(((n_A - 1) * var_A + (n_B - 1) * var_B) / (n_A + n_B - 2))
    cohen_d = (mu_A - mu_B) / (pooled_sd + 1e-16)
    return cohen_d


def get_DE_genes(
        gene_expression_matrix_df,
        cell_X_metadata_df,
        metadata_col,
        id0,
        id1,
        pseudocounts=1e-16,
        normalize_cell_counts=True,
        scale_factor=1e6,
        verbose: bool = False):
    '''
    To reproduce VISION tutorial FDR results, set normalize_cell_counts=True, normalization_factor=median counts,
    and pseudocounts=sqrt(n0 * n1) where n0 and n1 are the number of cells in each group. See:
    https://github.com/YosefLab/VISION/blob/master/R/Server.R#L704-L706
    Note that VISION does NOT normalize counts for you. You must do this yourself (as in the tutorial)
    Also, note that VISION uses Benjamini–Hochberg rather than Bonferroni to correct p-values, so
    you will get different p-values in VISION (and those aren't reproduced here - only FDR is reproduced).
    :param gene_expression_matrix_df: gene X cell df with counts.
    :param cell_X_metadata_df: cell X metadata df. index should be cell ids.
    :param metadata_col: column of cell_X_metadata_df to use for determining cell group
    :param id0: Value of metadata_col that identifies the first group
    :param id1: Value of metadata_col that identifies the second group
    :param pseudocounts: This number of pseudocounts will be added to each group during logFC computation
        E.g. if group 0 has mean 145 and group 1 has mean 0.03, then with 100 pseudocount, the logFC will
        change from log(145 / 0.034) ~ 8.5 to log((145 + 100) / (0.03 + 100)) ~ 0.9
    :param normalize_cell_counts: If True, each cell's counts will be normalized to sum to normalization_factor.
    :param scale_factor: Used only if normalize_cell_counts=True.
    '''
    if normalize_cell_counts:
        if gene_expression_matrix_df.isna().any().any():
            raise ValueError(f"gene_expression_matrix_df contains NaN's, I don't know how to normalize it!")
        gene_expression_matrix_df = \
            scale_factor * gene_expression_matrix_df.div(gene_expression_matrix_df.sum(axis=0), axis='columns')
    gene_names = list(gene_expression_matrix_df.index)
    scores = []
    avg_logFCs = []
    cohen_ds = []
    pct0s = []
    pct1s = []
    idx = cell_X_metadata_df.index  # Need this in case indices don't match! => .loc is your friend!
    group_0_all_cells = (cell_X_metadata_df.loc[idx, metadata_col] == id0)
    group_1_all_cells = (cell_X_metadata_df.loc[idx, metadata_col] == id1)
    group_0_or_1_all_cells = group_0_all_cells | group_1_all_cells
    for i, gene_name in enumerate(gene_names):
        # First subset the correct gene_counts and group indicators:
        # - Remove cells with nan counts
        # - Remove cells that are not group 0 or 1
        gene_counts = gene_expression_matrix_df.loc[gene_name, idx]
        is_not_nan = ~gene_counts.isna()
        use_cell = group_0_or_1_all_cells & is_not_nan
        gene_counts, group_0, group_1 = gene_counts[use_cell], group_0_all_cells[use_cell], group_1_all_cells[use_cell]
        # Done, now compute the test score!
        score = wilcoxon_score(gene_counts, group_0)
        scores.append(score)
        group_0_counts = gene_counts[group_0]
        group_1_counts = gene_counts[group_1]
        avg_logFCs.append(np.log((group_0_counts.mean() + pseudocounts) /
                                 (group_1_counts.mean() + pseudocounts)))
        cohen_ds.append(calc_cohens_d(group_0_counts, group_1_counts))
        pct0s.append(np.mean(group_0_counts > 0))
        pct1s.append(np.mean(group_1_counts > 0))
        if verbose and i % 1000 == 0:
            print(f"Processed {i} genes ... ")
    res_df = pd.DataFrame({'gene_name': gene_names,
                           'score': scores,
                           # 'p_val': p_vals, # This is computed at the end 'in parallel' for performance reasons.
                           'avg_logFC': avg_logFCs,
                           'cohen_d': cohen_ds,
                           'pct0': pct0s,
                           'pct1': pct1s,
                           # 'p_val_adj': p_vals_adj  # This is computed at the end 'in parallel' for performance.
                           })
    res_df['-|score|'] = -np.abs(res_df['score'])
    res_df['p_val'] = scipy.stats.norm.cdf(-np.abs(res_df['score'])) * 2.0
    res_df['p_val_adj'] = res_df['p_val'] * len(gene_names)
    res_df.sort_values(by='p_val_adj', inplace=True)
    res_df['rank'] = np.arange(1, res_df.shape[0] + 1, 1)
    res_df = res_df[["rank", "gene_name", "p_val", "avg_logFC", "cohen_d", "pct0",
                     "pct1", "p_val_adj", "score", "-|score|"]]
    return res_df
