r"""
The MeanModel imputes each feature using its mean value
"""
import numpy as np

from compass.turbo_mc.models.matrix_completion_model import MatrixCompletionModel


class MeanModel(MatrixCompletionModel):
    def __init__(self, n_epochs: int = 1):
        self.n_epochs = n_epochs
        self.column_means = None

    def _fit_matrix_init(
            self,
            X_observed: np.array,
            Z: None = None) -> None:
        self.nrows, self.ncols = X_observed.shape
        self.column_means = np.nanmean(X_observed, axis=0)

    def _step(self):
        pass

    def predict(self, r: int, c: int):
        r"""
        Predict entry (r, c)
        """
        return self.column_means[c]


class GlobalMeanModel(MatrixCompletionModel):
    def __init__(self, n_epochs: int = 1):
        self.n_epochs = n_epochs
        self.global_mean = None

    def _fit_matrix_init(
            self,
            X_observed: np.array,
            Z: None = None) -> None:
        self.nrows, self.ncols = X_observed.shape
        self.global_mean = np.nanmean(X_observed)

    def _step(self):
        pass

    def predict(self, r: int, c: int):
        r"""
        Predict entry (r, c)
        """
        return self.global_mean
