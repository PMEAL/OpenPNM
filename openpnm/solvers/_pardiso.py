from pypardiso import spsolve
from openpnm.solvers import DirectSolver
from scipy.sparse import csr_matrix, csc_matrix

__all__ = ['PardisoSpsolve']


class PardisoSpsolve(DirectSolver):
    """Solves a linear system using ``pypardiso.spsolve``."""

    def solve(self, A, b, **kwargs):
        """Solves the given linear system of equations Ax=b."""
        if not isinstance(A, (csr_matrix, csc_matrix)):
            A = A.tocsr()
        return (spsolve(A, b), 0)
