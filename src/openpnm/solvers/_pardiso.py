from scipy.sparse import csc_matrix, csr_matrix

from openpnm.solvers import DirectSolver

__all__ = ['PardisoSpsolve']


class PardisoSpsolve(DirectSolver):
    """Solves a linear system using ``pypardiso.spsolve``."""

    def solve(self, A, b, **kwargs):
        """Solves the given linear system of equations Ax=b."""
        from pypardiso import spsolve

        if not isinstance(A, (csr_matrix, csc_matrix)):
            A = A.tocsr()
        return (spsolve(A, b), 0)
