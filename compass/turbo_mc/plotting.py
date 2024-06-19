import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import Callable, Optional

from .matrix_manipulation import calc_is_nonconstant_column
from .metrics import compute_r2, compute_spearman_r2, compute_pearson_r2, compute_mse,\
    compute_spearman_r2s, compute_variance_explained, compute_variance_explained_fast
from .wilcoxon import get_DE_genes

pd.set_option('display.width', 1000)
pd.options.display.float_format = '{:,.2f}'.format


def plot_prediction(
        y_true: np.array,
        y_observed: np.array,
        y_completion: np.array,
        exclude_observed_for_plotting=True,
        exclude_observed_for_metric=True,
        metric='r2'):
    r"""
    Plot the true vs predicted value, annotated with the chosen metric.
    The metric shown excludes the observed entries.
    If the metric is mse, we give_each_column_equal_weight.
    """
    indices_to_use = np.where(np.isnan(y_observed)) if exclude_observed_for_plotting else\
        np.where(~np.isnan(y_true))  # Hacky way to get all entries...
    y_true_plotting, y_completion_plotting = \
        y_true[indices_to_use], y_completion[indices_to_use]
    plt.scatter(y_true_plotting, y_completion_plotting, s=8)
    minimum = min(np.min(y_true), np.min(y_completion)) - 1
    maximum = max(np.max(y_true), np.max(y_completion)) + 1
    plt.plot([minimum, maximum], [minimum, maximum], c='r')
    if metric == 'r2':
        r2 = compute_r2(y_true, y_observed, y_completion, exclude_observed=exclude_observed_for_metric)
        plt.title("r2 = %.2f" % r2)
    elif metric == 'mse':
        mse = compute_mse(y_true, y_observed, y_completion, exclude_observed=exclude_observed_for_metric,
                          give_each_column_equal_weight=True)
        plt.title("mse = %.2f" % mse)


def plot_predictions(
        X: np.array,
        X_observed: np.array,
        X_completion: np.array,
        k=5,
        figsize=(8, 6),
        exclude_observed=True,
        metric='r2',
        seed=1,
        exclude_observed_for_metric=True):
    r"""
    Plot k X k figure with random reactions predictions (true vs. predicted),
    annotated with the chosen metric. The metric shown excludes the observed entries.
    If the metric is mse, we give_each_column_equal_weight.
    If metric = 'reaction_idx', the plot will show the id (index) of the reaction in
    the plot, which is useful for e.g. following up on specific reactions for which
    the predictions look weird.
    """
    fig, a = plt.subplots(k, k, figsize=figsize)
    np.random.seed(seed)
    for i, r in enumerate(np.random.choice(range(X.shape[1]),
                          size=k * k, replace=False)):
        a[i // k][i % k].axis('off')
        indices_to_use = np.where(np.isnan(X_observed[:, r])) if exclude_observed else\
            np.where(~np.isnan(X[:, r]))  # Hacky way to get all entries...
        y_true = X[indices_to_use, r]
        y_pred = X_completion[indices_to_use, r]
        a[i // k][i % k].scatter(y_true, y_pred, s=8)
        minimum = min(np.min(y_true), np.min(y_pred)) - 1
        maximum = max(np.max(y_true), np.max(y_pred)) + 1
        a[i // k][i % k].plot([minimum, maximum], [minimum, maximum], c='r')
        if metric == 'r2':
            r2 = compute_r2(X[:, r], X_observed[:, r], X_completion[:, r],
                            exclude_observed=exclude_observed_for_metric)
            a[i // k][i % k].set_title("r2 = %.2f" % r2)
        elif metric == 'spearman_r2':
            spearman_r2 = compute_spearman_r2(X[:, r], X_observed[:, r], X_completion[:, r],
                                              exclude_observed=exclude_observed_for_metric)
            a[i // k][i % k].set_title("s_r2 = %.2f" % spearman_r2)
        elif metric == 'pearson_r2':
            pearson_r2 = compute_pearson_r2(X[:, r], X_observed[:, r], X_completion[:, r],
                                            exclude_observed=exclude_observed_for_metric)
            a[i // k][i % k].set_title("p_r2 = %.2f" % pearson_r2)
        elif metric == 'mse':
            mse = compute_mse(X[:, r], X_observed[:, r], X_completion[:, r],
                              exclude_observed=exclude_observed_for_metric,
                              give_each_column_equal_weight=True)
            a[i // k][i % k].set_title("mse = %.2f" % mse)
        elif metric == 'reaction_idx':
            a[i // k][i % k].set_title("r_idx = %d" % r)
    plt.show()


def compare_hypothesis_test_to_ground_truth(
        X_df_true: pd.DataFrame,
        X_df_pred: pd.DataFrame,
        hypothesis_test_true: pd.DataFrame,
        hypothesis_test_pred: pd.DataFrame,
        spearman_r2s_df: pd.DataFrame,
        spearman_r2s_mean: float,
        obs_metadata: pd.DataFrame,
        metadata_col: str,
        id0: str,
        id1: str,
        k: int = 5,
        figsize=(8, 12),
        seed: int = 1,
        show_largest_p_val_changes: bool = True,
        show_cohen_d: bool = True,
        show_top_features_in_any_direction: int = 20,
        show_top_features_in_each_directions: int = 0,
        feature_filter_func: Callable[[str], bool] = lambda x: True,
        cv_spearman_r2_df: Optional[pd.DataFrame] = None):
    r"""
    For pedanticness, X_df_true, X_df_pred, obs_metadata must only contain
    the observations from groups id0 and id1.
    :param X_df_true: True feature X obs DataFrame.
    :param X_df_pred: Predicted feature X obs DataFrame. MAY CONTAIN np.nan (e.g.
        for when performing the hypothesis test on only the observed entries.)
    :param hypothesis_test_true: True hypothesis test DataFrame; the output of get_DE_genes.
    :param hypothesis_test_pred: Predicted hypothesis test DataFrame; the output of get_DE_genes.
    :param spearman_r2s_df: DataFrame with spearman R2 for each feature.
        Used for showing the spearman R2 of each feature in the hypothesis test table,
        and also in the plots of the top extreme feature p-value changes. Although SR2s
        could be computed in this function, for performance reason they are passed to
        this method (not to recompute the spearman R2 of the ground truth multiple times
        when comparing against different hypothesis tests, e.g. MC completion HT and naive
        HT using observed entries). (NOTE: SR2s = Spearman R2, HT = Hypothesis Test).
    :param spearman_r2s_mean: Shown alongside the feature's spearman r2. It need not be
        the mean of spearman_r2s_df - indeed, it might e.g. exclude constant features.
    :param obs_metadata: The obs X metadata DataFrame.
    :param metadata_col: The name of the obs_metadata column to use
        for hypothesis testing (e.g. 'cell_type')
    :param id0: Value of metadata_col that identifies the first group
    :param id1: Value of metadata_col that identifies the second group
    :param k: A k * k grid of features will be plotted.
    :param figsize: Size of k * k figure
    :param seed: Numpy seed used to select k * k random features.
    :param show_largest_p_val_changes: If False, random features are plotted; else,
        the most extreme p-value changes will be plotted.
    :param show_cohen_d: If False, won't show cohen's d. (Which is useful to avoid
        wrapping around in print)
    :param show_top_features_in_any_direction: If provided, will show list of top
        differential features regardless of directionality. This can cause only
        features that point in one direction to dominate the list, whereas
        features in the other direction might be of interest.
    :param show_top_features_in_each_directions: If provided, will separately show
        top differential features that are greater/less.
    :param feature_filter_func: Function that, given a string (feature), determines if it
        should be included in the HT table. Of course, the feature can be excluded before calling
        compare_hypothesis_test_to_ground_truth, but sometimes it is better to filter deeper down.
    :param cv_spearman_r2_df: Cross-Validated Spearman R2 for each feature. If provided,
        will be shown in the HT table.
    """
    print(f"Merging hypothesis tests")
    hypothesis_tests =\
        pd.merge(
            hypothesis_test_true,
            hypothesis_test_pred,
            left_index=True,
            right_index=True
        )  # show_HT_in_one_direction operated on this object
    # Keep populating the hypothesis_tests object, then sanity check all data.
    hypothesis_tests['p_value_ratio'] =\
        np.log10(hypothesis_tests['p_val_x'] / hypothesis_tests['p_val_y'])
    hypothesis_tests['p_value_ratio_abs'] = np.abs(hypothesis_tests['p_value_ratio'])
    hypothesis_tests['-|score|_diff'] =\
        hypothesis_tests['-|score|_x'] - hypothesis_tests['-|score|_y']
    hypothesis_tests['-|score|_diff_abs'] = np.abs(hypothesis_tests['-|score|_diff'])

    def show_HT_in_one_direction(
            hypothesis_tests: pd.DataFrame,
            direction: Optional[int],
            rank_col: str,
            nfeatures: int):
        r"""
        :param hypothesis_tests: DataFrame which is merge of the groupd truth get_DE_genes
            and the predicted get_DE_genes. The merge must be such that the prefix '_x'
            corresponds to the true test, and '_y' corresponds to the predicted test.
            Must also have fields such as '-|score|_x', '-|score|_y', '-|score|_diff'
            populated.
        :param direction: If None, the top nfeatures features with lowest p-values will be
            shows. If 1, only the top nfeatures with 'cohen_d_x' > 0 will be shown. If -1,
            onle the top nfeatures with 'cohen_d_x' < 0 will be shown.
        :param rank_col: The column used to determine the ranking. Should be either 'rank_x' or
            'rank_y'. If it is 'rank_x', then the ground truth ranking is shown, else the
            predicted ranking is shown.
        :param nfeatures: How many top features to show.
        """
        if direction not in [1, -1, None]:
            raise ValueError(f"direction is '{direction}', not allowed.")
        if rank_col not in ['rank_x', 'rank_y']:
            raise ValueError(f"rank_col is '{rank_col}', not allowed.")
        true_or_imputed_str = 'True' if rank_col == 'rank_x' else 'Imputed'
        direction_str = 'mean(x) > mean(y)' if direction == 1 else 'mean(x) < mean(y)' if direction == -1 else ''
        print('*' * 10 + f" {true_or_imputed_str} top {nfeatures} features {direction_str} " + '*' * 10)
        # First pull out the features to use: either positive direction / negative direction / all.
        # Also filter based on custom function (e.g. for exclusion of exchange & transport reactions)
        features_to_use = (direction * hypothesis_tests['cohen_d_x']) > 0 if direction is not None\
            else [True] * hypothesis_tests.shape[0]
        features_to_use = features_to_use & np.array(list(map(feature_filter_func, hypothesis_tests.index)))
        sub_hypothesis_tests = hypothesis_tests[features_to_use].copy()
        # Create new ranks
        for what in ['rank_x', 'rank_y']:
            sub_hypothesis_tests.sort_values([what], inplace=True)
            sub_hypothesis_tests[what] = np.arange(1, sub_hypothesis_tests.shape[0] + 1, 1)
        # Report p-value deterioration metric
        print(f"%.2f pct p-values deteriorate (WARNING: Ties in ranks broken arbitrarily!)" %
              (sub_hypothesis_tests['-|score|_diff'] <= 0).mean())
        # Precisions at k
        top_true_tests = sub_hypothesis_tests.sort_values(['-|score|_x']).index
        top_imputed_tests = sub_hypothesis_tests.sort_values(['-|score|_y']).index
        for k0 in [1, 5, 10, 20, 50, 100]:
            precision_at_k0 = top_true_tests[:k0].isin(top_imputed_tests[:k0]).mean()
            print(f"precision at {k0} = {precision_at_k0}")
        print(sub_hypothesis_tests.sort_values([rank_col])[cols_to_use].head(nfeatures))

    # Sanity check that indices have same elements before adding Spearman R2 (plus other sanity checks)
    assert(X_df_true.shape == X_df_pred.shape)
    assert(X_df_true.shape[0] == hypothesis_test_true.shape[0])
    assert(X_df_pred.shape[0] == hypothesis_test_pred.shape[0])
    assert(hypothesis_tests.index.isin(spearman_r2s_df.index).all())
    assert(spearman_r2s_df.index.isin(hypothesis_tests.index).all())
    assert(hypothesis_tests.index.isin(X_df_true.index).all())
    assert(X_df_true.index.isin(hypothesis_tests.index).all())
    assert(hypothesis_tests.index.isin(X_df_pred.index).all())
    assert(X_df_pred.index.isin(hypothesis_tests.index).all())
    hypothesis_tests['spearman_r2'] = spearman_r2s_df['spearman_r2']  # Safe bc they have same index.
    if cv_spearman_r2_df is not None:
        assert(cv_spearman_r2_df.index.isin(spearman_r2s_df.index).all())
        assert(spearman_r2s_df.index.isin(cv_spearman_r2_df.index).all())
        hypothesis_tests['sr2_cv'] = cv_spearman_r2_df['cv_spearman_r2']  # Safe bc they have same index.

    # Print merged hypothesis tests (all features / positive features / negative features)
    cols_to_use = ['rank_x', 'rank_y', '-|score|_x', '-|score|_y', 'spearman_r2']
    if cv_spearman_r2_df is not None:
        cols_to_use += ['sr2_cv']
    if show_cohen_d:
        cols_to_use += ['cohen_d_x', 'cohen_d_y']
    if show_top_features_in_any_direction:  # All features
        nfeatures = show_top_features_in_any_direction
        # Print merged hypothesis test: top nfeatures
        show_HT_in_one_direction(hypothesis_tests, direction=None, rank_col='rank_x', nfeatures=nfeatures)
        show_HT_in_one_direction(hypothesis_tests, direction=None, rank_col='rank_y', nfeatures=nfeatures)
    if show_top_features_in_each_directions:  # positive and negative features
        nfeatures = show_top_features_in_each_directions
        show_HT_in_one_direction(hypothesis_tests, direction=1, rank_col='rank_x', nfeatures=nfeatures)
        show_HT_in_one_direction(hypothesis_tests, direction=1, rank_col='rank_y', nfeatures=nfeatures)
        show_HT_in_one_direction(hypothesis_tests, direction=-1, rank_col='rank_x', nfeatures=nfeatures)
        show_HT_in_one_direction(hypothesis_tests, direction=-1, rank_col='rank_y', nfeatures=nfeatures)

    print('*' * 10 + f" Global Stats " + '*' * 10)
    # Report p-value deterioration metric  # TODO: This is duplicated code
    print(f"%.2f pct p-values deteriorate" % (hypothesis_tests['-|score|_diff'] <= 0).mean())
    # Precisions at k  # TODO: This is duplicated code
    top_true_tests = hypothesis_tests.sort_values(['-|score|_x']).index
    top_imputed_tests = hypothesis_tests.sort_values(['-|score|_y']).index
    for k0 in [1, 5, 10, 20, 50, 100]:
        precision_at_k0 = top_true_tests[:k0].isin(top_imputed_tests[:k0]).mean()
        print(f"precision at {k0} = {precision_at_k0}")
    # Plot p-values scatterplot (currently I just plot the -|score| bc
    # some p-values are == 0.0)
    plt.title('True vs Imputed -|score|')
    x1, x2 = hypothesis_tests['-|score|_x'], hypothesis_tests['-|score|_y']
    plt.scatter(x1, x2, alpha=0.2)
    minimum = min(np.nanmin(x1), np.nanmin(x2)) - 1
    maximum = max(np.nanmax(x1), np.nanmax(x2)) + 1
    plt.plot([minimum, maximum], [minimum, maximum], c='r')
    plt.xlabel('True -|score|')
    plt.ylabel('Imputed -|score|')
    plt.show()
    # Plot most drastic p-value decreases upon completion (which can lead to false
    # positives)
    print(f"True vs imputed -|score| - Histogram of values per class:")
    fig, a = plt.subplots(2 * k, k, figsize=figsize)
    np.random.seed(seed)
    if show_largest_p_val_changes:
        least_deteriorated =\
            list(hypothesis_tests.sort_values(['-|score|_diff'], ascending=False).
                 index[: (k * k + 1) // 2])
        most_deteriorated =\
            list(hypothesis_tests.sort_values(['-|score|_diff'], ascending=False).
                 index[-(k * k - (k * k + 1) // 2):])
        generator = enumerate(least_deteriorated + most_deteriorated)
        print(f"=> Showing {len(least_deteriorated)} & {len(most_deteriorated)} most extreme "
              f"p-value decreases & increases after completion resp. (the former can lead to false positives!)")
    else:
        print(f"=> Showing {k * k} random features")
        generator = enumerate(np.random.choice(hypothesis_tests.index, size=k * k, replace=False))
    idx = obs_metadata.index
    group_0 = (obs_metadata.loc[idx, metadata_col] == id0)
    group_1 = (obs_metadata.loc[idx, metadata_col] == id1)
    if not (group_0 | group_1).all():
        raise ValueError(f"Not all observations belong to group 0 or 1. "
                         f"(The code can actually handle this case, but I prefer to be pedantic)")
    for i, feature in generator:
        vals_true = X_df_true.loc[feature, idx]
        vals_pred = X_df_pred.loc[feature, idx]
        a[2 * (i // k)][i % k].axis('off'), a[2 * (i // k) + 1][i % k].axis('off')
        a[2 * (i // k)][i % k].hist([vals_true[group_0], vals_true[group_1]])
        a[2 * (i // k)][i % k].set_title(
            "%s\n%.2f / %.2f" %
            (feature, hypothesis_tests.loc[feature, '-|score|_x'],
             hypothesis_tests.loc[feature, '-|score|_y']))
        # Subset the correct values and group indicators:
        # - Remove obs with nan vals
        is_not_nan = ~vals_pred.isna()
        vals_pred_notnan, group_0_notnan, group_1_notnan =\
            vals_pred[is_not_nan], group_0[is_not_nan], group_1[is_not_nan]
        a[2 * (i // k) + 1][i % k].hist([vals_pred_notnan[group_0_notnan], vals_pred_notnan[group_1_notnan]])
        a[2 * (i // k) + 1][i % k].set_title(
            "SR2: %.2f (%.2f)" %
            (spearman_r2s_df.loc[feature, 'spearman_r2'], spearman_r2s_mean))
    plt.tight_layout()
    plt.show()


def analyze_imputed_hypothesis_test(
        X_df_true: np.array,
        X_df_observed: np.array,
        X_df_pred: np.array,
        obs_metadata: pd.DataFrame,
        metadata_col: str,
        id0: str,
        id1: str,
        subset_obs: Optional[int] = 1000,
        seed: int = 1,
        show_cohen_d: bool = True,
        show_top_features_in_any_direction: int = 20,
        show_top_features_in_each_directions: int = 0,
        feature_filter_func: Callable[[str], bool] = lambda x: True,
        cv_spearman_r2_df: Optional[pd.DataFrame] = None,
        use_observed_in_completion: bool = False) -> None:
    r"""
    Given ground truth, observed, and predicted matrices, along with grouping
    metadata for each observation: performs the Wilcoxon rank-sum test based
    on the provided grouping, and reports tons of metrics and plots showing how
    the hypothesis testing results (e.g. top features with smallest p-values)
    changes when the observed and imputed matrices are used to conduct the
    hypothesis test instead of using the ground truth matrix.
    Only observations that are part of the two HT groups, and features of interest
    (as provided via feature_filter_func) are considered in this analysis. I.e.
    the first thing that happens is that we throw away all other observations and
    features.
    :param X_df_true: True feature X obs DataFrame.
    :param X_df_observed: Observed feature X obs DataFrame.
    :param X_df_pred: Imputed feature X obs DataFrame.
    :param obs_metadata: The obs X metadata DataFrame.
    :param metadata_col: The name of the obs_metadata column to use
        for hypothesis testing (e.g. 'cell_type')
    :param id0: Value of metadata_col that identifies the first group
    :param id1: Value of metadata_col that identifies the second group
    :param subset_obs: If provided, the X_df_true and X_df_pred hypothesis test
        will be performed on a subset of this number of observations.
        This greatly speeds up test computation when there are many observations.
    :param seed: Seed for all sources of randomness.
    :param show_cohen_d: If False, won't show cohen's d. (Which is useful to avoid
        wrapping around in print)
    :param show_top_features_in_any_direction: If provided, will show list of top
        differential features regardless of directionality. This can cause only
        features that point in one direction to dominate the list, whereas
        features in the other direction might be of interest.
    :param show_top_features_in_each_directions: If provided, will separately show
        top differential features that are greater/less.
    :param feature_filter_func: Function that, given a string (feature), determines if it
        should be included in the HT table. Of course, the feature can be excluded before calling
        analyze_imputed_hypothesis_test, but sometimes it is better to filter deeper down so that
        e.g. variance explained is not affected by upstream filtering.
    :param cv_spearman_r2_df: Cross-Validated Spearman R2 for each feature. If provided,
        will be shown in the HT table.
    :param use_observed_in_completion: If True, will use the observed values in the completion,
        thus overriding whatever the model predicted. (Some models already do this).
    """
    # Subset the observations from each group first, so the correct spearman R2 shows.
    print(f"Original matrix size: {X_df_true.shape}")
    print(f"Original percent of matrix observed: {(~np.isnan(X_df_observed.to_numpy())).mean()} / 1.00")
    print(f"Subsetting only observations from the two groups involved in the hypothesis test: "
          f"'{metadata_col}' in ['{id0}', '{id1}'], and only features of interest.")
    features_to_use = np.array(list(map(feature_filter_func, X_df_true.index)))
    obs_to_subset = obs_metadata.index[obs_metadata[metadata_col].isin([id0, id1])]
    X_df_true = X_df_true.loc[features_to_use, obs_to_subset]
    X_df_observed = X_df_observed.loc[features_to_use, obs_to_subset]
    X_df_pred = X_df_pred.loc[features_to_use, obs_to_subset]
    if use_observed_in_completion:
        indices_to_replace = ~(X_df_observed.isna())
        X_df_pred[indices_to_replace] = X_df_observed[indices_to_replace]
    if cv_spearman_r2_df is not None:
        cv_spearman_r2_df = cv_spearman_r2_df.loc[features_to_use]
    obs_metadata = obs_metadata.loc[obs_to_subset]
    print(f"New matrix size: {X_df_true.shape}")
    print(f"New percent of matrix observed: {(~np.isnan(X_df_observed.to_numpy())).mean()} / 1.00")
    # Recover the obs X features matrices.
    X = X_df_true.transpose().to_numpy()
    X_observed = X_df_observed.transpose().to_numpy()
    X_completion = X_df_pred.transpose().to_numpy()
    # Report variance explained
    print("Variance explained fast (pseudo_std=1e-1) = %.2f"
          % compute_variance_explained_fast(X, X_completion, pseudo_std=1e-1))
    print("Variance explained (pedantic) = %.2f" % compute_variance_explained(X, X_completion))
    # Plot predictions for some features
    print('*' * 10 + f" Plotting completions for 25 random features " + '*' * 10)
    plot_predictions(
        X=X,
        X_observed=X,
        X_completion=X_completion,
        exclude_observed=False,
        seed=seed,
        exclude_observed_for_metric=False,
        metric='spearman_r2')
    # Plot spearman_r2 histogram, for non-constant features
    is_nonconstant_column = calc_is_nonconstant_column(X)
    spearman_r2s = compute_spearman_r2s(X, X, X_completion, exclude_observed=False)
    plt.title('Histogram of Spearman R2s per feature\nexcluding constant features.')
    plt.hist(spearman_r2s[is_nonconstant_column], bins=30)
    plt.show()
    spearman_r2s_mean = spearman_r2s[is_nonconstant_column].mean()
    spearman_r2s_df = pd.DataFrame({'feature': X_df_true.index, 'spearman_r2': spearman_r2s})
    spearman_r2s_df.set_index(['feature'], inplace=True)
    print(f"Mean Spearman R2 (excluding constant columns) = %.2f" % spearman_r2s_mean)
    # Plot histogram of number of observations per column
    plt.title(f"Histogram of observations per column")
    plt.hist((~np.isnan(X_observed)).sum(axis=0), bins=20)
    plt.show()

    def hypothesis_test_func(X: np.array, subset_obs: Optional[int] = None):
        if subset_obs:
            print(f"WARNING: Performing hypothesis test on a subset of "
                  f"{min(X.shape[1], subset_obs)} observations.")
            np.random.seed(seed)
            random_obs =\
                np.random.choice(
                    obs_metadata.index,
                    min(X.shape[1], subset_obs),
                    replace=False)
        else:
            random_obs = obs_metadata.index
        res = get_DE_genes(
            gene_expression_matrix_df=X.loc[:, random_obs],
            cell_X_metadata_df=obs_metadata.loc[random_obs],
            metadata_col=metadata_col,
            id0=id0,
            id1=id1,
            normalize_cell_counts=False)
        res.rename(columns={'gene_name': 'feature'}, inplace=True)
        res.set_index(['feature'], inplace=True)
        return res
    print(f"Perfoming HT on true data ...")
    hypothesis_test_true = hypothesis_test_func(X_df_true, subset_obs)
    print(f"Perfoming HT on observed data ...")
    hypothesis_test_observed = hypothesis_test_func(X_df_observed, None)
    print(f"Perfoming HT on imputed data ...")
    hypothesis_test_pred = hypothesis_test_func(X_df_pred, subset_obs)

    def compare_hypothesis_test_to_ground_truth_func(
            X_df_pred,
            hypothesis_test_pred,
            spearman_r2s_df,
            spearman_r2s_mean):
        return compare_hypothesis_test_to_ground_truth(
                X_df_true,
                X_df_pred,
                hypothesis_test_true,
                hypothesis_test_pred,
                spearman_r2s_df,
                spearman_r2s_mean,
                obs_metadata,
                metadata_col,
                id0,
                id1,
                seed=seed,
                show_cohen_d=show_cohen_d,
                show_top_features_in_any_direction=show_top_features_in_any_direction,
                show_top_features_in_each_directions=show_top_features_in_each_directions,
                feature_filter_func=feature_filter_func,
                cv_spearman_r2_df=cv_spearman_r2_df)

    print('*' * 50)
    print('*' * 15 + ' COMPLETION METRICS ' + '*' * 15)
    print('*' * 50)
    compare_hypothesis_test_to_ground_truth_func(X_df_pred, hypothesis_test_pred,
                                                 spearman_r2s_df, spearman_r2s_mean)
    print('*' * 51)
    print('*' * 15 + ' OBSERVATION METRICS ' + '*' * 15)
    print('*' * 51)
    compare_hypothesis_test_to_ground_truth_func(X_df_observed, hypothesis_test_observed,
                                                 spearman_r2s_df * 0 + 1.0, 1.0)
