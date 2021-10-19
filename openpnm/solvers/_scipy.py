from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.linalg import spsolve
from openpnm.solvers import DirectSolver


class ScipySpsolve(DirectSolver):

    def solve(self, A, b, **kwargs):
        if not isinstance(A, (csr_matrix, csc_matrix)):
            A = A.tocsr()
        return (spsolve(A, b), 0)
