import sys
import time
from abc import ABC, abstractmethod
import copy
from dataclasses import dataclass
import logging
import numpy as np
from typing import List, Optional, Tuple, Union, Callable

from compass.turbo_mc.iterative_models.matrix_oracle import MatrixOracle
from compass.turbo_mc.models.matrix_completion_model import MatrixCompletionModel
from compass.turbo_mc.models.cv_matrix_completion_model import CVMatrixCompletionModel
from compass.turbo_mc.matrix_manipulation import get_list_of_random_matrix_indices

from compass import globals


class IterativeMatrixCompletionModel(ABC):  # pragma: no cover
    @abstractmethod
    def fit(self, matrix_oracle: MatrixOracle, Z: Optional[np.array] = None):
        raise NotImplementedError

    @abstractmethod
    def predict_all(self) -> np.array:
        raise NotImplementedError

    @abstractmethod
    def observed_matrix(self) -> np.array:
        raise NotImplementedError

    @abstractmethod
    def cv_spearman_r2s(self) -> np.array:
        raise NotImplementedError

    def impute_all(self) -> np.array:
        r"""
        Impute all unobserved matrix entries (not inplace).
        This is a handy top-class method.
        """
        X_observed, X_predicted = self.observed_matrix(), self.predict_all()
        unobserved_indices = np.where(np.isnan(X_observed))
        X_imputed = X_observed.copy()
        X_imputed[unobserved_indices] = X_predicted[unobserved_indices]
        return X_imputed


class IterativeMCMWithPrescribedPcts(IterativeMatrixCompletionModel):
    r"""
    Given the percent of observations to perform in each round
    (sampling_densities), does iterative querying.
    """
    def __init__(
            self,
            sampling_densities: List[float],
            cv_models: List[CVMatrixCompletionModel],
            verbose: bool = False):
        if len(sampling_densities) != len(cv_models):
            raise ValueError('sampling_densities and cv_models should have the same length!')
        self.sampling_densities = sampling_densities
        self.cv_models = copy.deepcopy(cv_models)
        self.verbose = verbose

    def fit(self, matrix_oracle: MatrixOracle, Z: Optional[np.array] = None):
        X_observed = matrix_oracle.observed_matrix()
        if not np.all(np.isnan(X_observed)):
            raise ValueError(
                "Received a warm-started oracle, but I don't yet have logic to handle "
                "that case smartly.")
        curr_cv_spearman_r2s = None  # type: Optional[np.array]
        for i, (sampling_density, cv_model) in enumerate(zip(self.sampling_densities, self.cv_models)):
            if self.verbose:
                print(f"IterativeMCMWithPrescribedPcts: Fitting iteration {i}: sampling_density = {sampling_density}")
            if curr_cv_spearman_r2s is None:
                # First round of fitting, just sample uniformly
                R, C = matrix_oracle.shape()
                matrix_indices = get_list_of_random_matrix_indices(R, C, sampling_density)
            else:
                # We have CV information, use it!
                matrix_indices = smart_list_of_matrix_indices(R, C, X_observed, sampling_density, curr_cv_spearman_r2s)
            matrix_oracle.observe_entries(matrix_indices)
            X_observed = matrix_oracle.observed_matrix()
            cv_model.fit_matrix(X_observed, Z)
            curr_cv_spearman_r2s = cv_model.cv_spearman_r2s()
        self.X_observed = matrix_oracle.observed_matrix()
        self.X_completion = cv_model.predict_all()
        self.curr_cv_spearman_r2s = curr_cv_spearman_r2s
        if self.verbose:
            sampling_density = (~np.isnan(self.X_observed)).mean()
            print(f"IterativeMCMWithPrescribedPcts: Finished fitting. Final sampling density: {sampling_density}")

    def predict_all(self) -> np.array:
        return self.X_completion

    def cv_spearman_r2s(self) -> np.array:
        return self.curr_cv_spearman_r2s

    def observed_matrix(self) -> np.array:
        return self.X_observed


def smart_list_of_matrix_indices(
        R: int,
        C: int,
        X_observed: np.array,
        sampling_density: float,
        cv_spearman_r2s: np.array
) -> List[Tuple[int, int]]:
    r"""
    Queries the bottom half of features with poorest spearman R2.
    """
    assert(len(cv_spearman_r2s) == C)
    columns_with_poor_performance = np.argsort(cv_spearman_r2s)[: (C // 2)]
    observations_per_column = int(2 * R * sampling_density)
    res = []  # type: List[Tuple[int, int]]
    for c in columns_with_poor_performance:
        unobserved_rows = np.where(np.isnan(X_observed[:, c]))[0]
        chosen_rows = np.random.choice(unobserved_rows, observations_per_column, replace=False)
        res += [(r, c) for r in chosen_rows]
    return res


@dataclass
class IterativeMCMState:
    R: int
    C: int
    sampled_density: float


class IterativeMCMWithGuaranteedSpearmanR2(IterativeMatrixCompletionModel):
    r"""
    Iteratively queries the matrix until _every single column_ has
    cross-validated spearman R2 above a requested threshold.
    Hard-to-predict columns might be fully queried as a result.

    This guarantees that all columns have high-quality reconstructions!

    The predicted matrix (i.e. the completion) uses the observed
    entry if possible, else imputes the entry using the underlying
    low-rank model.
    """
    def __init__(
            self,
            cv_model: Union[CVMatrixCompletionModel,
                            Callable[[IterativeMCMState], CVMatrixCompletionModel]],
            requested_cv_spearman_r2: float = 0.7,
            sampling_density: float = 0.01,
            finally_refit_model: Optional[Union[MatrixCompletionModel,
                                                Callable[[IterativeMCMState], MatrixCompletionModel]]] = None,
            min_pct_meet_sr2_requirement: float = 1.0,
            verbose: bool = False,
            plot_progress: bool = False,
            max_iterations: int = 1000000,
            logger_dir: str = ""):
        r"""
        :param cv_model: What CV model is fit at each step. If a different model
            wants to be fit at each iteration, a function mapping iteration number
            to a CVMatrixCompletionModel can be provided instead.
        :param requested_cv_spearman_r2: The requested Spearman R2. Every column
            of the imputed matrix is guaranteed (in terms of CV performance) to
            have this generalization Spearman R2.
        :param sampling_density: What percent of the matrix is observed at each step.
            I.e., the budget we have at each iteration.
        :param finally_refit_model: If provided, a final unique model is fit on
            all observations (This typically works better than the cv_model's
            ensemble predictions). Since we don't know beforehand at what iteration
            of the process we will stop, a function that takes the iteration number and
            returns a MatrixCompletionModel can be provided instead.
        :param min_pct_meet_sr2_requirement: How many of the columns should meet the
            requested_cv_spearman_r2. Set to 1.0 means _all_ columns must meet.
        :param plot_progress: If True, will plot histogram of CV Spearman after each round.
        :param max_iterations: By default the model will query the matrix until the minimum
            spearman R2 condition is met. If you want to make sure that no more than certain
            number of iterations are made, you can set this argument.
        """
        self.cv_model_func = cv_model if callable(cv_model) else lambda state: copy.deepcopy(cv_model)
        self.requested_cv_spearman_r2 = requested_cv_spearman_r2
        self.sampling_density = sampling_density
        self.finally_refit_model_func = finally_refit_model
        self.min_pct_meet_sr2_requirement = min_pct_meet_sr2_requirement
        self.verbose = verbose
        self.plot_progress = plot_progress
        self.max_iterations = max_iterations
        globals.init_logger(logger_dir)
        self.logger = logging.getLogger("compass")

    def fit(self, matrix_oracle: MatrixOracle, Z: Optional[np.array] = None):
        fit_start_time = time.time()
        # R: # of cells
        # C: # of reactions
        R, C = matrix_oracle.shape()
        sampling_density = self.sampling_density
        X_observed = matrix_oracle.observed_matrix()
        curr_cv_spearman_r2s = None  # type: Optional[np.array]
        for iteration in range(self.max_iterations):
            iteration_start_time = time.time()
            self.logger.info('\n' + '*' * 8 + f" Iteration {iteration + 1} " + '*' * 8)
            cv_model =\
                self.cv_model_func(
                    IterativeMCMState(
                        R=R,
                        C=C,
                        sampled_density=sampling_density * (iteration + 1)))

            # for each reaction, choose which cells to observe entry
            if iteration == 0:
                # First round of fitting!
                if np.all(np.isnan(X_observed)):
                    # Matrix is completely unobserved: just sample uniformly
                    matrix_indices = get_list_of_random_matrix_indices(R, C, sampling_density)
                else:
                    # The oracle is already warm-started: we want to start by knowing what those entries are!
                    observed_rows, observed_cols = np.where(~np.isnan(X_observed))
                    matrix_indices = list(zip(observed_rows, observed_cols))
            else:
                # We have CV information, use it!
                matrix_indices =\
                    _choose_entries_from_underperforming_columns(
                        X_observed,
                        sampling_density,
                        curr_cv_spearman_r2s,
                        self.requested_cv_spearman_r2)
            self.logger.info("Querying MatrixOracle ...")
            matrix_oracle.observe_entries(matrix_indices, iteration)
            X_observed = matrix_oracle.observed_matrix()
            # Fit cv_model
            pct_observed_so_far = (~np.isnan(X_observed)).mean()
            self.logger.info("Percent of matrix observed so far: %.3f/1.00" % pct_observed_so_far)
            self.logger.info("Fitting CVMatrixCompletionModel ...")
            cv_model.fit_matrix(X_observed, Z)
            # Recompute CV Spearman R2s
            curr_cv_spearman_r2s = cv_model.cv_spearman_r2s()
            _override_fully_observed_columns_sr2_to_1(curr_cv_spearman_r2s, X_observed)
            worse_cv_spearman_r2 = np.sort(curr_cv_spearman_r2s)[
                min(int(C * (1.0 - self.min_pct_meet_sr2_requirement)), C - 1)]
            self.logger.info(
                "Worse@%.2f vs Mean vs Requested CV Spearman R2 = %.3f vs. %.3f vs %.3f"
                % (self.min_pct_meet_sr2_requirement, worse_cv_spearman_r2, curr_cv_spearman_r2s.mean(),
                   self.requested_cv_spearman_r2))
            if self.plot_progress:  # pragma: no cover
                plt.title("Histogram of current CV Spearman R2s")
                plt.hist(curr_cv_spearman_r2s, bins=20)
                plt.show()
                plt.title("Histogram of observations per column")
                plt.hist((~np.isnan(X_observed)).sum(axis=0), bins=20)
                plt.show()
            self.logger.info("Iteration finished! Time: %.0fs" % (time.time() - iteration_start_time))
            if worse_cv_spearman_r2 >= self.requested_cv_spearman_r2:
                if self.verbose:
                    self.logger.info(
                        "Achieved goal of %.3f CV Spearman R2 >= %.3f!"
                        % (self.min_pct_meet_sr2_requirement, self.requested_cv_spearman_r2))
                break
        # The reported cv_spearman_r2s are the ones of the last CV model.
        self.curr_cv_spearman_r2s = curr_cv_spearman_r2s
        X_observed = matrix_oracle.observed_matrix()
        self.X_observed = X_observed
        # Refit one final large model if requested.
        if self.finally_refit_model_func is not None:
            if callable(self.finally_refit_model_func):
                self.logger.info(
                    f"\nRefitting final model on ({R}, {C}) matrix after observing "
                    + "%.3f/1.00 pct of entries ..." % (sampling_density * (iteration + 1)))
                final_model =\
                    self.finally_refit_model_func(
                        IterativeMCMState(
                            R=R,
                            C=C,
                            sampled_density=sampling_density * (iteration + 1)))
            else:
                self.logger.info("\nNot refitting any model at the very end.")
                final_model = copy.deepcopy(self.finally_refit_model_func)
            final_model.fit_matrix(X_observed, Z)
        else:
            final_model = cv_model
        self.X_predicted = final_model.predict_all()
        pct_observed = (~np.isnan(X_observed)).mean()
        self.logger.info(
            f"Finished iterative fitting! "
            f"Number of iterations used: {iteration + 1}; "
            + "Pct of matrix sampled: %.3f/1.00; " % pct_observed
            + "Total time: %.0fs" % (time.time() - fit_start_time))

    def predict_all(self) -> np.array:
        r"""
        You'll want to use impute_all instead, since some columns might have
        been fully observed because they were hard to predict.
        """
        return self.X_predicted

    def cv_spearman_r2s(self) -> np.array:
        return self.curr_cv_spearman_r2s

    def observed_matrix(self) -> np.array:
        return self.X_observed


def _override_fully_observed_columns_sr2_to_1(
        curr_cv_spearman_r2s: np.array,
        X_observed: np.array) -> None:
    r"""
    The curr_cv_spearman_r2s of fully observed columns is made 1.
    Note that the np.array curr_cv_spearman_r2s is modified in-place.
    """
    is_fully_observed = (~np.isnan(X_observed)).all(axis=0)
    curr_cv_spearman_r2s[is_fully_observed] = 1.0


def _choose_entries_from_underperforming_columns(
        X_observed: np.array,
        sampling_density: float,
        cv_spearman_r2s: np.array,
        requested_cv_spearman_r2: float
) -> List[Tuple[int, int]]:
    r"""
    Distributes budget evenly amoung all columns whose spearman R2 is below the
    requested threshold. Any remaining budget is distributed among the remaining
    columns.
    :param X_observed: Currently observed matrix
    :param sampling_density: How many entries we seek to sample.
    :param cv_spearman_r2s: The cross-validated Spearman R2 of each column
    :param requested_cv_spearman_r2: The minimum Spaerman R2 we want to achieve.
    :return: The list of (row, column) pairs where to sample next based on this strategy.
    """
    assert(len(cv_spearman_r2s[cv_spearman_r2s <= requested_cv_spearman_r2]) > 0)
    R, C = X_observed.shape
    assert(len(cv_spearman_r2s) == C)
    unobserved_per_column = np.isnan(X_observed).sum(axis=0)
    total_budget = int(sampling_density * R * C)
    # print(f"total_budget = {total_budget}")
    res = []  # type: List[Tuple[int, int]]
    unobserved_per_poor_column = unobserved_per_column[cv_spearman_r2s <= requested_cv_spearman_r2]
    n_unobserved_for_poor_columns = unobserved_per_poor_column.sum()
    # print(f"n_unobserved_for_poor_columns = {n_unobserved_for_poor_columns}")
    if n_unobserved_for_poor_columns >= total_budget:
        # print(f"Can allocate all resources amongst poor columns.")
        # Binary search for n_obs_per_column
        # print(f"Binary searching for number of observations per column.")
        lo = 0  # Too little
        hi = R  # Enough
        while lo < hi - 1:
            mid = (lo + hi) // 2
            if np.minimum(unobserved_per_poor_column, mid).sum() < total_budget:
                lo = mid
            else:
                hi = mid
        n_obs_per_column = hi
        # print(f"n_obs_per_column = {n_obs_per_column}")
        columns_sorted_by_sr2 = np.argsort(cv_spearman_r2s)
        for c in columns_sorted_by_sr2:
            # Select entries for this column uniformly at random.
            unobserved_rows = np.where(np.isnan(X_observed[:, c]))[0]
            n_obs_for_this_column = min(min(n_obs_per_column, len(unobserved_rows)), total_budget)
            chosen_rows = np.random.choice(unobserved_rows, n_obs_for_this_column, replace=False)
            total_budget -= n_obs_for_this_column
            res += [(r, c) for r in chosen_rows]
        assert(total_budget == 0)
        # print(f"res = {res}")
    else:
        # print(f"Will query all poor columns, then distribute remaining queries amongst all other columns")
        column_still_available = [True] * C
        for c in range(C):
            if cv_spearman_r2s[c] <= requested_cv_spearman_r2:
                for r in range(R):
                    if np.isnan(X_observed[r, c]):
                        res.append((r, c))
                column_still_available[c] = False
        assert(len(res) == n_unobserved_for_poor_columns)
        total_budget = total_budget - len(res)
        # Split remaining budget uniformly amongst 'good' columns.
        rows, cols = np.where(np.isnan(X_observed))
        remaining_entries = list(zip(rows, cols))
        remaining_entries = [(r, c) for (r, c) in remaining_entries if column_still_available[c]]
        # print(f"remaining_entries = {remaining_entries}")
        # if len(remaining_entries) <= total_budget:
        #     print(f"I ended up querying the full matrix!")
        indices = np.random.choice(len(remaining_entries),
                                   size=min(len(remaining_entries), total_budget),
                                   replace=False)
        # print(f"indices = {indices}")
        res += [remaining_entries[i] for i in indices]
        # print(f"res = {res}")
    return res
