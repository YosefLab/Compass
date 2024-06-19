r"""
This module contains the interface CVMatrixCompletionModel and its concretions.
A CVMatrixCompletionModel is a MatrixCompletionModel which can additionally
respond to the 'cv_spearman_r2s' method.

The concretion KFoldCVMatrixCompletionModel wrap a given
MatrixCompletionModel. Said model is fit once for each fold. Out-of fold
predictions are computed and used to compute out-of-fold Spearman R2 for
each feature (column).
"""
import sys
import time
from abc import ABC, abstractmethod
import copy
import logging
import numpy as np
import scipy.stats
from typing import Optional

from compass.turbo_mc.metrics import are_equal, is_constant
from compass.turbo_mc.models.matrix_completion_model import MatrixCompletionModel

# Logger for KFoldCVMatrixCompletionModel
logger = logging.getLogger(__name__ + ".KFoldCVMatrixCompletionModel")
logger.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s %(lineno)s %(name)s: %(message)s")
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(formatter)
logger.addHandler(handler)
logger = None

# Logger for TrainValSplitCVMatrixCompletionModel
logger = logging.getLogger(__name__ + ".TrainValSplitCVMatrixCompletionModel")
logger.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s %(lineno)s %(name)s: %(message)s")
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(formatter)
logger.addHandler(handler)
logger = None


def compute_spearman_r2(
        y_true: np.array,
        y_pred: np.array) -> float:
    r"""
    Standard Spearman R2 correlation: a number between -1 and 1.
    Note that:
        - predicted columns that are almost equal to the ground truth are declared Spearman R2 1.0
        - predicted columns that are not almost equal yet ground truth is constant, have Spearman R2 0.0.
        - for all other columns, the usual Spearman R2 is computed
    """
    assert(len(y_true) == len(y_pred))
    # Border case 0: there are no true labels or predictions. Return 0.0 (random)
    if len(y_true) == 0:
        return 0.0
    # Border case 1: predictions are very close to ground truth
    if are_equal(y_true, y_pred):
        return 1.0
    # Border case 2: predictions are not very close to ground truth, yet ground truth is constant
    if not are_equal(y_true, y_pred) and is_constant(y_true):
        return 0.0
    # Usual computation:
    # Border case: all true values or predicted values are constant (but the other is not)
    if np.std(y_true) == 0.0:  # pragma: no cover
        raise ValueError("This should never happen!")
    if np.std(y_pred) == 0.0:
        assert(np.std(y_true) != 0.0)
        return 0.0
    spearman_r2 = scipy.stats.spearmanr(y_true, y_pred).correlation
    # Even then, sometimes spearman_r2 can return np.nan :( Not sure why
    if np.isnan(spearman_r2):  # pragma: no cover
        # print("WARNING: scipy.stats.spearmanr returned np.nan")
        return 0.0
    return spearman_r2


def compute_spearman_r2s(
        X_true: np.array,
        X_pred: np.array)\
        -> np.array:
    r"""
    Computes Spearman R2 accross all columns.
    Note that:
        - predicted columns that are almost equal to the ground truth are declared Spearman R2 1.0
        - predicted columns that are not almost equal yet ground truth is constant, have Spearman R2 0.0.
        - for all other columns, the usual Spearman R2 is computed
    """
    def validate_input():
        assert(X_true.shape == X_pred.shape)
        true_nan_indices = np.where(np.isnan(X_true))
        pred_nan_indices = np.where(np.isnan(X_pred))
        assert((true_nan_indices[0] == pred_nan_indices[0]).all())
        assert((true_nan_indices[1] == pred_nan_indices[1]).all())
    validate_input()
    nrows, ncols = X_true.shape
    spearman_r2s = np.zeros(shape=(ncols))
    for c in range(ncols):
        not_nan_rows = np.where(~np.isnan(X_true[:, c]))
        y_true = X_true[not_nan_rows, c].flatten()
        y_pred = X_pred[not_nan_rows, c].flatten()
        spearman_r2 = compute_spearman_r2(y_true, y_pred)
        spearman_r2s[c] = spearman_r2
    return spearman_r2s


class CVMatrixCompletionModel(MatrixCompletionModel, ABC):
    @abstractmethod
    def cv_spearman_r2s(self):  # pragma: no cover
        raise NotImplementedError


class KFoldCVMatrixCompletionModel(CVMatrixCompletionModel):
    def __init__(
            self,
            model: MatrixCompletionModel,
            n_folds: int = 5,
            verbose: bool = False,
            finally_refit_model: Optional[MatrixCompletionModel] = None):
        r"""
        :param finally_refit_model: If provided, the model is fit on the full data
            at the very end. This tends to produce a better fit than the
            ensemble of smaller models. Useful for a final model fit.
            Usless if you only care about the cross-valudation metrics.
        """
        self.model = copy.deepcopy(model)
        self.n_folds = n_folds
        self.verbose = verbose
        self.finally_refit_model = copy.deepcopy(finally_refit_model)
        self.logger = logging.getLogger(__name__ + ".KFoldCVMatrixCompletionModel")

    def fit_matrix(self, X_observed: np.array, Z: Optional[np.array] = None):
        start_time = time.time()
        X_completion = np.zeros_like(X_observed)
        X_completion_cv = np.zeros_like(X_observed) + np.nan
        n_folds = self.n_folds
        # TODO: My CV-split is simple and efficient, but might (with very low prob) assign all
        # observed entries in a column to the same fold, which I think would blow up downstream code.
        # An improvement would be to make sure the observations in each column are equally distributed
        # among all folds instead.
        X_fold_numbers = np.random.randint(low=0, high=n_folds, size=X_observed.shape)
        X_fold_numbers[np.where(np.isnan(X_observed))] = -1
        for fold_number in range(n_folds):
            self.logger.info(f"Fitting fold {fold_number + 1} ...")
            model_fold = copy.deepcopy(self.model)
            X_observed_fold = np.copy(X_observed)
            X_observed_fold[np.where(X_fold_numbers == fold_number)] = np.nan
            model_fold.fit_matrix(X_observed_fold, Z)
            X_completion_fold = model_fold.predict_all()
            X_completion += X_completion_fold / n_folds
            out_of_fold_indices = np.where(X_fold_numbers == fold_number)
            X_completion_cv[out_of_fold_indices] = X_completion_fold[out_of_fold_indices]
        if self.finally_refit_model:
            self.logger.info("Refitting model on all data ...")
            full_model = copy.deepcopy(self.finally_refit_model)
            full_model.fit_matrix(X_observed, Z)
            X_completion = full_model.predict_all()
        # Save variables
        self.X_completion = X_completion
        self.X_completion_cv = X_completion_cv
        self.X_fold_numbers = X_fold_numbers
        self.X_observed = X_observed
        self.logger.info("Finished fitting CV model! Total time: %.0fs" % (time.time() - start_time))

    def predict(self, r: int, c: int) -> float:
        return self.X_completion[r, c]

    def predict_all(self) -> np.array:
        return self.X_completion

    def cv_spearman_r2s(self) -> np.array:
        return compute_spearman_r2s(self.X_observed, self.X_completion_cv)


class DummyCVMatrixCompletionModel(CVMatrixCompletionModel):  # pragma: no cover
    r"""
    The CV Spearman R2 is simply the training Spearman R2! (i.e. far from unbiased!)
    """
    def __init__(
            self,
            model: MatrixCompletionModel):
        self.model = copy.deepcopy(model)

    def fit_matrix(self, X_observed: np.array, Z: Optional[np.array] = None):
        self.model.fit_matrix(X_observed, Z)
        self.X_completion = self.model.predict_all()
        self.X_completion_cv = np.copy(self.X_completion)
        self.X_completion_cv[np.where(np.isnan(X_observed))] = np.nan
        self.X_observed = X_observed

    def predict(self, r: int, c: int) -> float:
        return self.X_completion[r, c]

    def predict_all(self) -> np.array:
        return self.X_completion

    def cv_spearman_r2s(self) -> np.array:
        return compute_spearman_r2s(self.X_observed, self.X_completion_cv)


class TrainValSplitCVMatrixCompletionModel(CVMatrixCompletionModel):
    def __init__(
            self,
            model: MatrixCompletionModel,
            train_ratio: float = 0.8,
            verbose: bool = False,
            finally_refit_model: Optional[MatrixCompletionModel] = None):
        r"""
        :param train_ratio: What percentage of the data to assign to training
            (the rest is assigned to the validation set).
        :param finally_refit_model: If provided, the model is fit on the full data
            at the very end. This tends to produce a better fit than the
            ensemble of smaller models. Useful for a final model fit.
            Usless if you only care about the cross-valudation metrics.
        """
        self.model = copy.deepcopy(model)
        self.train_ratio = train_ratio
        self.verbose = verbose
        self.finally_refit_model = copy.deepcopy(finally_refit_model)
        self.logger = logging.getLogger(__name__ + ".TrainValSplitCVMatrixCompletionModel")

    def fit_matrix(self, X_observed: np.array, Z: Optional[np.array] = None):
        start_time = time.time()
        train_ratio = self.train_ratio
        model = self.model
        R, C = X_observed.shape
        is_train = np.zeros(shape=(R, C)) != 0  # I.e. all False
        for c in range(C):
            observed_row_indices = np.where(~np.isnan(X_observed[:, c]))[0]
            train_indices_for_this_col =\
                np.random.choice(
                    observed_row_indices,
                    int(len(observed_row_indices) * train_ratio),
                    replace=False)
            is_train[train_indices_for_this_col, c] = True
        is_valid = ~is_train
        # Remove all unobserved entries from train and valid.
        is_train = is_train & ~np.isnan(X_observed)
        is_valid = is_valid & ~np.isnan(X_observed)
        # Create training data.
        X_train = X_observed.copy() + np.nan
        X_train[np.where(is_train)] = X_observed[np.where(is_train)]
        # Fit model to training data.
        self.logger.info("Training model ...")
        model.fit_matrix(X_train, Z)
        # Make predictions
        X_prediction = model.predict_all()
        # Make predictions for validation set only.
        X_valid_pred = X_prediction + np.nan
        X_valid_pred[is_valid] = X_prediction[is_valid]
        # Make targets for validation set only
        X_valid = X_observed.copy() + np.nan
        X_valid[is_valid] = X_observed[is_valid]
        # If finally_refit_model, do so, and upate the predictions.
        if self.finally_refit_model:
            self.logger.info("Refitting model on all data ...")
            full_model = copy.deepcopy(self.finally_refit_model)
            full_model.fit_matrix(X_observed, Z)
            X_prediction = full_model.predict_all()
        # Save results
        self.X_prediction = X_prediction
        self.X_valid = X_valid
        self.X_valid_pred = X_valid_pred
        self.logger.info("Finished fitting CV model! Total time: %.0fs" % (time.time() - start_time))

    def predict(self, r: int, c: int) -> float:
        return self.X_prediction[r, c]

    def predict_all(self) -> np.array:
        return self.X_prediction

    def cv_spearman_r2s(self) -> np.array:
        return compute_spearman_r2s(self.X_valid, self.X_valid_pred)
