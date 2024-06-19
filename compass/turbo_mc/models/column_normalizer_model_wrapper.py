import numpy as np
import time
from typing import List

from compass.turbo_mc.models.matrix_completion_model import MatrixCompletionModel
from compass.turbo_mc.utils import timeit  # For profiling
import compass.turbo_mc.utils as utils  # For profiling


class ColumnNormalizerModelWrapper(MatrixCompletionModel):
    r"""
    Normalizes the columns of the matrix before fitting the underlying MatrixCompletionModel.
    """
    def __init__(self, model: MatrixCompletionModel):
        self.n_epochs = model.n_epochs
        self.model = model

    def _fit_matrix_init(
            self,
            X_observed: np.array,
            Z: None = None) -> None:
        if utils.TIMEIT:  # For profiling
            print("\t\tColumnNormalizerModelWrapper::_fit_matrix_init ...")  # For profiling
            time_start = time.time()  # For profiling
        self.column_stds = np.nanstd(X_observed, axis=0)
        self.column_stds[np.isnan(self.column_stds)] = 1.0  # Deals with empty columns
        self.column_means = np.nanmean(X_observed, axis=0)
        self.column_means[np.isnan(self.column_means)] = 0.0  # Deals with empty columns
        X_observed_normalized = (X_observed - self.column_means) / (self.column_stds + 1e-16)
        if utils.TIMEIT:  # For profiling
            time_tot = time.time() - time_start  # For profiling
            print("\t\tColumnNormalizerModelWrapper::_fit_matrix_init, time = %.2f" % time_tot)  # For profiling
        self.model._fit_matrix_init(X_observed_normalized, Z)

    def _step(self):
        self.model._step()

    def objective(self, regularized: bool = True) -> float:
        return self.model.objective(regularized)

    def objective_upper_bounds(self) -> List[float]:
        return self.model.objective_upper_bounds()

    def predict(self, r: int, c: int) -> float:
        r"""
        Predict entry (r, c).
        Implemented by querying the underlying model, then un-normalizing the result.
        :param r: row
        :param c: column
        :return: prediction
        """
        return self.model.predict(r, c) * self.column_stds[c] + self.column_means[c]

    @timeit("\t\tColumnNormalizerModelWrapper")
    def _after_model_predict_call(self, res):
        return res * np.expand_dims(self.column_stds, axis=0) + np.expand_dims(self.column_means, axis=0)

    def predict_all(self) -> np.array:
        r"""
        Predict all matrix entries (even those observed).
        """
        res = self.model.predict_all()
        res = self._after_model_predict_call(res)
        return res
