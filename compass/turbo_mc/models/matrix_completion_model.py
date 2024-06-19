r"""
Abstract base class for all Matrix Completion Models. A Matrix Completion Model
must implement the following two methods:
- fit
- predict
Additionally, the predict_all method should be overriden when possible with a
faster implementation, since the default 'predict_all' implementation just loops
over all matrix entries and calls 'predict'.
"""
from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from typing import Optional

from compass.turbo_mc.matrix_manipulation import matrix_from_observations, get_observations_from_partially_observed_matrix


class MatrixCompletionModel(ABC):
    def fit(
            self,
            observations: pd.DataFrame,
            nrows: int,
            ncols: int,
            Z: Optional[np.array] = None) -> None:
        r"""
        Fits the model to the observed data, given optional additional row features Z.
        (Z conveying 'latent space')
        :param observations: pd.DataFrame of triples (row, col, val)
        :param nrows: Number of rows of the matrix.
        :param ncols: Number of columns of the matrix.
        :param Z: 2D np.array of optional additional row features.
        """
        self._fit_init(observations, nrows, ncols, Z)
        for epoch in range(self.n_epochs):
            self._step()

    def fit_matrix(self, X_observed: np.array, Z: Optional[np.array] = None) -> None:
        r"""
        Fits the model to the observed data, given optional additional row features Z.
        (Z conveying 'latent space')
        :param X_observed: 2D matrix of observations X features.
        :param Z: 2D np.array of optional additional row features.
        """
        self._fit_matrix_init(X_observed, Z)
        for epoch in range(self.n_epochs):
            self._step()

    def _fit_init(
            self,
            observations: pd.DataFrame,
            nrows: int,
            ncols: int,
            Z: Optional[np.array] = None) -> None:
        X_observed = matrix_from_observations(observations, nrows=nrows, ncols=ncols)
        self._fit_matrix_init(X_observed, Z)

    def _fit_matrix_init(
            self,
            X_observed: np.array,
            Z: None = None) -> None:
        observations = get_observations_from_partially_observed_matrix(X_observed)
        self._fit_init(observations, X_observed.shape[0], X_observed.shape[1], Z)

    @abstractmethod
    def predict(self, r: int, c: int) -> float:
        r"""
        Predict a specific entry of the matrix.
        :param r: row
        :param c: column
        :return: prediction
        """
        raise NotImplementedError  # pragma: no cover

    def predict_all(self) -> np.array:
        r"""
        Predict _all_ entries of the unobserved matrix using the underlying model.
        Observed matrix entries should still be computed from the model, _not_ memorized.
        If you just want to impute the _missing_ entries use the impute_all method instead.
        NOTE: To use this default top-class method, the underlying model must have nrows and ncols set.
        TODO: Wrap the 'fit' method to auto-populate nrows and ncols.
        """
        if self.nrows is None or self.ncols is None:
            raise KeyError("self.nrows or self.ncols not initialized, cannot predict_all")  # pragma: no cover
        res = np.zeros(shape=(self.nrows, self.ncols))
        for r in range(self.nrows):
            for c in range(self.ncols):
                res[r, c] = self.predict(r, c)
        return res

    def impute_all(self, X_observed: np.array, verbose: bool = False) -> np.array:
        r"""
        Impute all unobserved matrix entries (not inplace).
        This is a handy top-class method.
        :param X_observed: Observed matrix. Unobserved entries have np.nan.
        """
        if verbose:
            print("In impute_all ... ")  # pragma: no cover
        unobserved_indices = np.where(np.isnan(X_observed))
        X_completion = X_observed.copy()
        X_completion[unobserved_indices] = self.predict_all()[unobserved_indices]
        if verbose:
            print("Out impute_all")  # pragma: no cover
        return X_completion
