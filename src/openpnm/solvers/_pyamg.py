import pyamg
from scipy.sparse import csr_matrix

from ._base import IterativeSolver

__all__ = ['PyamgRugeStubenSolver']


# write a PyamgRugeStubenSolver class
class PyamgRugeStubenSolver(IterativeSolver):
    """Iterative solver based on PyAMG's `ruge_stuben_solver`."""

    def solve(self, A, b, x0=None):
        if not isinstance(A, csr_matrix):
            A = A.tocsr()
        ml = pyamg.ruge_stuben_solver(A)
        return ml.solve(b, x0=x0, tol=self.tol, maxiter=self.maxiter, return_info=True)
