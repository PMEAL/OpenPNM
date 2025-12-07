from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.linalg import spsolve, cg
from openpnm.solvers import DirectSolver, IterativeSolver

__all__ = ['ScipySpsolve', 'ScipyCG']


class ScipySpsolve(DirectSolver):
    """Solves a linear system using ``scipy.sparse.linalg.spsolve``."""

    def solve(self, A, b, **kwargs):
        """Solves the given linear system of equations Ax=b."""
        if not isinstance(A, (csr_matrix, csc_matrix)):
            A = A.tocsr()
        return (spsolve(A, b), 0)


class ScipyCG(IterativeSolver):
    """Solves a linear system using ``scipy.sparse.linalg.cg``."""

    def solve(self, A, b, **kwargs):
        """Solves the given linear system of equations Ax=b."""
        if not isinstance(A, (csr_matrix, csc_matrix)):
            A = A.tocsr()
        atol = self._get_atol(b)
        try:
            return cg(A, b, tol=self.tol, atol=atol, **kwargs)
        except TypeError:
            return cg(A, b, rtol=self.tol, atol=atol, **kwargs)
