import numpy as np
import time
from typing import Dict, List

from compass.turbo_mc.models.matrix_completion_model import MatrixCompletionModel
from compass.turbo_mc.matrix_manipulation import get_nonconstant_columns
from compass.turbo_mc.utils import timeit  # for profiling
import compass.turbo_mc.utils as utils  # for profiling


class ExcludeConstantColumnsModelWrapper(MatrixCompletionModel):
    r"""
    Excludes constant columns when fitting the underlying model.
    Constant columns are predicted with their constant value.
    This wrapper is usefull for two reasons:
    1) Fitting the underlying model is faster as the training data is smaller.
    2) Constant columns are guaranteed to be predicted as constant, which fixes the r2 for those columns.
        (Fitting models such as MF may lead to some noise in constant columns, and as such, -1.0 r2s)
    """
    def __init__(self, model: MatrixCompletionModel):
        self.n_epochs = model.n_epochs
        self.model = model

    def _fit_matrix_init(
            self,
            X_observed: np.array,
            Z: None = None) -> None:
        if utils.TIMEIT:  # For profiling
            print("\tExcludeConstantColumnsModelWrapper::_fit_matrix_init ...")  # For profiling
            time_start = time.time()  # For profiling
        self.nrows, self.ncols = X_observed.shape
        self.column_means = np.nanmean(X_observed, axis=0)
        self.column_means[np.isnan(self.column_means)] = 0.0  # Deals with empty columns
        self.nonconstant_column_indices = get_nonconstant_columns(X_observed)
        self.is_nonconstant_column_index = [False] * self.ncols
        for idx in self.nonconstant_column_indices:
            self.is_nonconstant_column_index[idx] = True
        self.map_from_original_index_to_new_index = dict()  # type: Dict[int, int]
        for i in range(len(self.nonconstant_column_indices)):
            idx_of_nonconstant_column = self.nonconstant_column_indices[i]
            self.map_from_original_index_to_new_index[idx_of_nonconstant_column] = i
        if len(self.nonconstant_column_indices) == 0:
            self.all_columns_are_constant = True
        else:
            self.all_columns_are_constant = False
            X_observed_subset = X_observed[:, self.nonconstant_column_indices]
            if utils.TIMEIT:  # For profiling
                time_tot = time.time() - time_start  # For profiling
                print("\tExcludeConstantColumnsModelWrapper::_fit_matrix_init, time = %.2f"
                      % time_tot)  # For profiling
            self.model._fit_matrix_init(X_observed_subset, Z=Z)

    def _step(self):
        if not self.all_columns_are_constant:
            self.model._step()

    def objective(self, regularized: bool = True) -> float:
        return self.model.objective(regularized)

    def objective_upper_bounds(self) -> List[float]:
        return self.model.objective_upper_bounds()

    def predict(self, r: int, c: int) -> float:
        r"""
        If the column is constant, returns that constant value, else queries the
        underlying model.
        :param r: row
        :param c: column
        :return: prediction
        """
        if not self.is_nonconstant_column_index[c]:
            return self.column_means[c]
        else:
            new_c = self.map_from_original_index_to_new_index[c]
            return self.model.predict(r, new_c)

    @timeit("\tExcludeConstantColumnsModelWrapper")
    def _after_model_predict_call(self, predictions_for_nonconstant_columns):
        res = np.ones(shape=(self.nrows, 1)) @ np.expand_dims(self.column_means, axis=0)
        res[:, self.nonconstant_column_indices] = predictions_for_nonconstant_columns
        return res

    def predict_all(self) -> np.array:
        r"""
        Predict all matrix entries (even those observed).
        """
        if not self.all_columns_are_constant:
            predictions_for_nonconstant_columns = self.model.predict_all()
            res = self._after_model_predict_call(predictions_for_nonconstant_columns)
        else:
            res = np.ones(shape=(self.nrows, 1)) @ np.expand_dims(self.column_means, axis=0)
        return res
