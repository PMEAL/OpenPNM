from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.linalg import spsolve
from openpnm.solvers import DirectSolver

__all__ = ['ScipySpsolve']


class ScipySpsolve(DirectSolver):
    """Brief description of 'ScipySpsolve'"""

    def solve(self, A, b, **kwargs):
        """Brief description of 'solve'"""
        if not isinstance(A, (csr_matrix, csc_matrix)):
            A = A.tocsr()
        return (spsolve(A, b), 0)
