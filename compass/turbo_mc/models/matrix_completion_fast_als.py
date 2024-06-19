r"""
Matrix Completion based on maximum-margin matrix factorization (MMMF), termed softImputeALS:
It is Algorithm 5.1 in:
`Matrix Completion and Low-Rank SVD via Fast Alternating Least Squares`
http://jmlr.org/papers/volume16/hastie15a/hastie15a.pdf

Solves the non-convex optimization problem:
min_{A, B} 1/2 * |P_Omega(X) - P_Omega(AB^T)|_F^2 + lam/2 * |A|_F^2 + lam/2 * |B|_F^2

This is essentially just ALS but densifying the matrix before solving the LS problems,
which leads to them having the same structure and so can be solved way faster
(vectorized) rather than looping over all rows/columns.

Note that we do not implement the sparse version (since in my intended use case the matrices
are not so sparse), nor do we implement the orthogonalization and redistribution in step
3 of Algorithm 3.1, nor step 5 of Algorithm 3.1. These two simplifications to Algorithm 3.1
mean that the original objective is not guaranteed to decrease at each iteration, but is not
so relevant in practice so we exclude them.
"""
import numpy as np

from compass.turbo_mc.models.matrix_completion_model import MatrixCompletionModel
from compass.turbo_mc.utils import timeit


class MatrixCompletionFastALS(MatrixCompletionModel):
    def __init__(
            self,
            n_factors: int,
            lam: float,
            n_epochs: int,
            verbose: bool):
        self.n_factors = n_factors
        self.lam = lam
        self.n_epochs = n_epochs
        self.verbose = verbose

    def _zero_out(self, X, indices):
        res = X.copy()
        res[indices] = 0
        return res

    @timeit("\t\t\tMatrixCompletionFastALS")
    def _fit_matrix_init(
            self,
            X_observed: np.array,
            Z: None = None) -> None:
        nrows, ncols = X_observed.shape[0], X_observed.shape[1]
        observed_indices = np.where(~np.isnan(X_observed))
        unobserved_indices = np.where(np.isnan(X_observed))
        A = np.random.normal(size=(nrows, self.n_factors)) / np.sqrt(self.n_factors)
        B = np.random.normal(size=(ncols, self.n_factors)) / np.sqrt(self.n_factors)
        X_zeroed_out = self._zero_out(X_observed, unobserved_indices)
        # Cache variables
        self.unobserved_indices = unobserved_indices
        self.X_observed, self.observed_indices, self.A, self.B, self.X_zeroed_out =\
            X_observed, observed_indices, A, B, X_zeroed_out
        self.epoch = 0

    @timeit("\t\t\tMatrixCompletionFastALS")
    def _step(self):
        if self.verbose:
            print(f"Epoch {self.epoch}")
        # Recover variables
        observed_indices, A, B, X_zeroed_out, zero_out =\
            self.observed_indices, self.A, self.B, self.X_zeroed_out, self._zero_out
        # A update
        X_star = X_zeroed_out + zero_out(A @ B.T, observed_indices)
        A = X_star @ B @ np.linalg.inv(B.T @ B + self.lam * np.eye(B.shape[1]))
        # B update
        X_star = X_zeroed_out + zero_out(A @ B.T, observed_indices)
        B = X_star.T @ A @ np.linalg.inv(A.T @ A + self.lam * np.eye(A.shape[1]))
        # Compute matrix completion
        self.X_completion = A @ B.T
        # Cache variables
        self.A, self.B = A, B
        self.epoch += 1

    @timeit("\t\t\tMatrixCompletionFastALS")
    def objective(self, regularized: bool = True) -> float:
        r"""
        Returns the value of the optimization objective:
        1/2 * |P_Omega(X) - P_Omega(AB^T)|_F^2 + lam/2 * |A|_F^2 + lam/2 * |B|_F^2
        :param regularized: If False, only 1/2 * |P_Omega(X) - P_Omega(Z)|_F^2 is returned.
        """
        X_observed, X_completion, unobserved_indices, A, B, lam, zero_out =\
            self.X_observed, self.X_completion, self.unobserved_indices, self.A, self.B, self.lam, self._zero_out
        res = 0.5 * np.linalg.norm(zero_out(X_observed - X_completion, unobserved_indices), 'fro') ** 2
        if regularized:
            res += 0.5 * lam * np.linalg.norm(A, 'fro') ** 2 + 0.5 * lam * np.linalg.norm(B, 'fro') ** 2
        return res

    def predict(self, r: int, c: int) -> float:
        return self.X_completion[r, c]

    @timeit("\t\t\tMatrixCompletionFastALS")
    def predict_all(self) -> np.array:
        r"""
        Predict ALL the matrix (including the observed entries). Returns matrix.
        Convenient for computing the full MSE between the ground truth and the low rank
        approximation - this MSE is lower bounded by the truncated SVD's MSE.
        """
        return self.X_completion.copy()
