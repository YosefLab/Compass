r"""
'MatrixOracle' is the abstract base class for all matrix oracles. A
matric oracle provides a means of observing an underlying matrix by
means of an 'observe_entries' method, as well as a 'shape' method
that gives the dimension of the matrix.

Subclasses only have to implement 'shape' and the helper method '_observe_entries'.
"""
import numpy as np
from typing import List, Tuple, Iterable

from abc import ABC, abstractmethod


class MatrixOracle(ABC):
    @abstractmethod
    def shape(self) -> Tuple[int, int]:  # pragma: no cover
        r"""
        Returns the shape of the underlying matrix.
        """
        raise NotImplementedError

    def observe_entries(
        self,
        matrix_indices: List[Tuple[int, int]],
        iter: int
    ) -> np.array:
        if 'X_observed' not in self.__dict__:  # TODO: This is duplicated code
                self.X_observed = np.zeros(shape=self.shape()) + np.nan
        # Extract rows and cols; make them np.arrays for easier indexing.
        rows, cols = map(np.array, zip(*matrix_indices))
        # Exclude (row, col) pairs that were already observed (no need to re-compute them!)
        index_not_observed = np.isnan(self.X_observed[(rows, cols)])
        unobserved_rows, unobserved_cols = rows[index_not_observed], cols[index_not_observed]
        # Query the unobserved entries (cast rows and cols back to List[int] to respect interface)
        unobserved_vals = self._observe_entries(list(unobserved_rows), list(unobserved_cols), iter)
        self.X_observed[(unobserved_rows, unobserved_cols)] = unobserved_vals
        return self.X_observed[(rows, cols)]

    @abstractmethod
    def _observe_entries(
        self,
        rows: List[int],
        cols: List[int],
        iter: int
    ) -> np.array:   # pragma: no cover
        r"""
        Returns the one-dimensional np.array of values found at the requested entries
        (as specified by 'rows' and 'cols').
        """
        raise NotImplementedError

    def observed_matrix(self) -> np.array:
        r"""
        Returns the currently observed matrix.
        """
        if 'X_observed' not in self.__dict__:  # TODO: This is duplicated code
            self.X_observed = np.zeros(shape=self.shape()) + np.nan
        return self.X_observed

    def _add_observations(
        self,
        rows: List[int],
        cols: List[int],
        vals: List[float]
    ) -> None:
        r"""
        Method that subclasses can call to add observations to the matrix.
        Handy for e.g. warm-starting an oracle upon __init__().
        """
        assert(len(rows) == len(cols))
        assert(len(rows) == len(vals))
        if 'X_observed' not in self.__dict__:  # TODO: This is duplicated code
            self.X_observed = np.zeros(shape=self.shape()) + np.nan
        self.X_observed[(rows, cols)] = vals


class OracleWithAPrescribedMatrix(MatrixOracle):
    def __init__(self, X: np.array):
        self.X = X

    def shape(self) -> Tuple[int, int]:
        return self.X.shape

    def _observe_entries(
        self,
        rows: List[int],
        cols: List[int]
    ) -> np.array:
        return self.X[(rows, cols)]


class WarmstartedOracleWithAPrescribedMatrix(MatrixOracle):
    def __init__(
        self,
        X: np.array,
        warm_started_indices: List[Tuple[int, int]]
    ):
        self.X = X
        rows, cols = zip(*warm_started_indices)
        vals = X[(rows, cols)]
        self._add_observations(rows, cols, vals)

    def shape(self) -> Tuple[int, int]:
        return self.X.shape

    def _observe_entries(
        self,
        rows: List[int],
        cols: List[int]
    ) -> np.array:
        return self.X[(rows, cols)]
