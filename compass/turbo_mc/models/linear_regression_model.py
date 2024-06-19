r"""
The LinearRegressionModel regresses each matrix column against the row features in Z.
This is the same as saying we have one linear regression problem per matrix column.
"""
import numpy as np
import pandas as pd
import sklearn.linear_model

from .matrix_completion_model import MatrixCompletionModel


class LinearRegressionModel(MatrixCompletionModel):
    def __init__(
            self,
            alpha: float = 0.0,
            fit_intercept: bool = False,
            solver: str = 'sklearn',
            standardize_features: bool = False):
        r"""
        :param alpha: L2 regularization strength
        :param fit_intercept: If to fit intercept to model.
        :param solver: The solver to used, by default 'sklearn', which uses
            sklearn.linear_model.Ridge.
        :param standardize_features: If True, features will be standardized to
            mean 0 and std 1 before fitting the lienar regression model.
        """
        self.n_epochs = 1
        self.alpha = alpha
        self.fit_intercept = fit_intercept
        known_solvers = ['sklearn', 'manual']
        if solver not in known_solvers:
            raise ValueError(f"Unknown solver: {solver}. Available solvers: {known_solvers}")
        if solver == 'manual' and fit_intercept:
            raise ValueError('Manual solver does not support fitting the intercept.')
        self.solver = solver
        self.standardize_features = standardize_features

    def _fit_init(
            self,
            observations: pd.DataFrame,
            nrows: int,
            ncols: int,
            Z: np.array) -> None:
        self.nrows = nrows
        self.ncols = ncols
        self.nfactors = Z.shape[1]
        self.Z = Z.copy()
        if self.standardize_features:
            self.Z = (self.Z - self.Z.mean(axis=0, keepdims=True)) / (self.Z.std(axis=0, keepdims=True) + 1e-16)
        self.column_embeddings = np.zeros(shape=(self.nfactors, self.ncols))
        self.intercepts = np.zeros(shape=(1, self.ncols))
        self.observations_by_column = observations.groupby('col')

    def _step(self):
        Z, observations_by_column = self.Z, self.observations_by_column
        for c in range(self.ncols):
            # Solve the L2-regularized linear regression problem for column c
            rows_observed_for_this_column = list(observations_by_column.get_group(c)['row'])

            # Build linear system Aw = b
            A = Z[rows_observed_for_this_column, :]
            b = np.array(list(observations_by_column.get_group(c)['val']))

            # Solve system for w (and intercept)
            if self.solver == 'sklearn':
                # Solve with sklearn
                model = \
                    sklearn.linear_model.Ridge(
                        alpha=self.alpha,
                        fit_intercept=self.fit_intercept)\
                    .fit(A, b)
                w = model.coef_
                intercept = model.intercept_
            elif self.solver == 'manual':
                w = np.linalg.inv(A.T @ A + self.alpha * np.eye(A.shape[1])) @ A.T @ np.reshape(b, newshape=(-1, 1))
                w = w.flatten()
                intercept = 0
            else:
                raise ValueError(f"Unknown solver {self.solver}")  # pragma: no cover

            self.column_embeddings[:, c] = w
            self.intercepts[0, c] = intercept

    def predict(self, r: int, c: int) -> float:
        r"""
        Predict a specific entry of the matrix.
        :param r: row
        :param c: column
        :return: prediction
        """
        return self.Z[r, :] @ self.column_embeddings[:, c] + self.intercepts[0, c]

    def predict_all(self) -> np.array:
        r"""
        Predict all entries of the matrix using the underlying linear model.
        """
        return self.Z @ self.column_embeddings + self.intercepts
