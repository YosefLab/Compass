r"""
Matrix Factorization Model wrapping SVD from the surprise library.

Optimizes the objective:
1/2 * |P_Omega(X) - P_Omega(AB^T)|_F^2 + lam_A/2 * |A|_F^2 + 0.5 * lam_B/2 * |B|_F^2
where lam_A = reg_all * |Omega| / R and lam_B = reg_all * |Omega| / C

Allows for row/column biased with the 'biased=True' flag upon initialization.
"""
import numpy as np
import pandas as pd
import surprise

from compass.turbo_mc.models.matrix_completion_model import MatrixCompletionModel
from compass.turbo_mc.matrix_manipulation import matrix_from_observations


class MatrixFactorizationModel(MatrixCompletionModel):
    def __init__(
            self,
            n_factors,
            random_state=42,
            reg_all=0.0,
            n_epochs=3000,
            lr_all=0.03,
            biased=False,
            verbose=False):
        self.n_epochs = n_epochs
        self.reg_all = reg_all
        self.model = surprise.SVD(
            n_factors=n_factors,
            random_state=random_state,
            reg_all=reg_all,
            n_epochs=n_epochs,
            lr_all=lr_all,
            biased=biased,
            verbose=False)  # I'll take care of the verbose
        self.verbose = verbose

    def _zero_out(self, X, indices):
        res = X.copy()
        res[indices] = 0
        return res

    def _fit_init(
            self,
            observations: pd.DataFrame,
            nrows: int,
            ncols: int,
            Z: None = None):
        self.nrows = nrows
        self.ncols = ncols
        reader = surprise.Reader(rating_scale=(-1e16, 1e16))
        observations_for_surprise = observations.rename(columns={
            "row": "userID",
            "col": "itemID",
            "val": "rating"
        })
        data = surprise.Dataset.load_from_df(observations_for_surprise, reader)
        trainset = data.build_full_trainset()
        self.trainset = trainset
        self.model._fit_init(trainset)
        self.epoch = 0
        X_observed = matrix_from_observations(observations, nrows, ncols)
        self.X_observed = X_observed
        self.unobserved_indices = np.where(np.isnan(X_observed))
        self.omega_size = observations.shape[0]

    def _step(self):
        if self.verbose:
            print(f"Epoch {self.epoch}")
        self.model._step()
        # Compute completion
        self.X_completion = self.get_row_embeddings() @ self.get_col_embeddings()
        if self.model.biased:
            self.X_completion += self.get_row_biases() + self.get_col_biases() + self.get_global_mean()
        self.epoch += 1

    def objective(self, regularized: bool = True) -> float:
        r"""
        Returns the value of the optimization objective:
        1/2 * |P_Omega(X) - P_Omega(AB^T)|_F^2 + lam_A/2 * |A|_F^2 + 0.5 * lam_B/2 * |B|_F^2
        (If the model is biased, AB^T must be ajusted with the biases)
        :param regularized: If False, only 1/2 * |P_Omega(X) - P_Omega(AB^T)|_F^2 is returned.
        """
        X_observed, X_completion, unobserved_indices, reg_all, omega_size, A, B =\
            self.X_observed, self.X_completion, self.unobserved_indices, self.reg_all, self.omega_size, self.model.pu,\
            self.model.qi
        res = 0.5 * np.linalg.norm(self._zero_out(X_observed - X_completion, unobserved_indices), 'fro') ** 2
        if regularized:
            lam_A = reg_all * omega_size / X_observed.shape[0]
            lam_B = reg_all * omega_size / X_observed.shape[1]
            res += 0.5 * lam_A * np.linalg.norm(A, 'fro') ** 2 + 0.5 * lam_B * np.linalg.norm(B, 'fro') ** 2
        return res

    def predict(self, r: int, c: int) -> float:
        return self.model.predict(r, c).est

    def get_row_embeddings(self) -> np.array:
        r"""
        returns the matrix of row embeddings, of size nrows X n_factors
        """
        nrows = self.model.pu.shape[0]
        row_permutation = [self.trainset.to_inner_uid(i) for i in range(nrows)]
        return self.model.pu[row_permutation, :]

    def get_row_biases(self) -> np.array:
        nrows = self.model.pu.shape[0]
        row_permutation = [self.trainset.to_inner_uid(i) for i in range(nrows)]
        res = np.expand_dims(self.model.bu[row_permutation], axis=1)
        return res

    def get_col_embeddings(self) -> np.array:
        r"""
        returns the matrix of column embeddings, of size ncols X n_factors
        """
        ncols = self.model.qi.shape[0]
        column_permutation = [self.trainset.to_inner_iid(i) for i in range(ncols)]
        return np.transpose(self.model.qi[column_permutation, :])

    def get_col_biases(self) -> np.array:
        ncols = self.model.qi.shape[0]
        column_permutation = [self.trainset.to_inner_iid(i) for i in range(ncols)]
        res = np.expand_dims(self.model.bi[column_permutation], axis=0)
        return res

    def get_global_mean(self) -> float:
        return self.trainset.global_mean

    def predict_all(self) -> np.array:
        r"""
        Predict ALL the matrix (including the observed entries). Returns matrix.
        Convenient for computing the full MSE between the ground truth and the low rank
        approximation - this MSE is lower bounded by the truncated SVD's MSE.
        """
        return self.X_completion.copy()
