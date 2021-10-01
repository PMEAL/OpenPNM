import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.linalg import spsolve
from openpnm.solvers import DirectSolver, Integrator


class ScipySpsolve(DirectSolver):

    def solve(self, A, b, **kwargs):
        if not isinstance(A, (csr_matrix, csc_matrix)):
            A = A.tocsr()
        return (spsolve(A, b), 0)


class ScipyRK45(Integrator):

    def __init__(self, atol=1e-6, rtol=1e-6, verbose=False, linsolver=None):
        self.atol = atol
        self.rtol = rtol
        self.verbose = verbose
        self.linsolver = linsolver

    def solve(self, rhs, x0, tspan, saveat, **kwargs):
        options = {
            "atol": self.atol,
            "rtol": self.rtol,
            "t_eval": saveat,
            "verbose": self.verbose,
        }
        sol = solve_ivp(rhs, tspan, x0, method="RK45", **options)
        if sol.success:
            return TransientSolution(sol.t, sol.y)
        raise Exception(sol.message)


class Solution(np.ndarray): ...


class TransientSolution(Solution):

    def __new__(cls, t, x):
        obj = np.asarray(x).view(cls)
        obj.t = np.asarray(t).view(cls)
        return obj

    def _create_interpolant(self):
        self._interpolant = interp1d(self.t, self, bounds_error=True)

    def interpolate(self, t):
        if not hasattr(self, "_interpolant"):
            self._create_interpolant()
        return self._interpolant(t)

    __call__ = interpolate
