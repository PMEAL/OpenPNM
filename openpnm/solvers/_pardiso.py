from pypardiso import spsolve
from openpnm.solvers import DirectSolver
from scipy.sparse import csr_matrix, csc_matrix


class PardisoSpsolve(DirectSolver):

    def solve(self, A, b, **kwargs):
        if not isinstance(A, (csr_matrix, csc_matrix)):
            A = A.tocsr()
        # TODO: solver.solve should return (x, info) not just x
        return (spsolve(A, b), 0)
