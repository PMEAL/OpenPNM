from pypardiso import spsolve
from openpnm.solvers import DirectSolver
from scipy.sparse import csr_matrix, csc_matrix

__all__ = ['PardisoSpsolve']


class PardisoSpsolve(DirectSolver):
    """Brief description of 'PardisoSpsolve'"""

    def solve(self, A, b, **kwargs):
        """Brief description of 'solve'"""
        if not isinstance(A, (csr_matrix, csc_matrix)):
            A = A.tocsr()
        # TODO: solver.solve should return (x, info) not just x
        return (spsolve(A, b), 0)
