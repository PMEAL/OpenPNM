from scipy.integrate import solve_ivp

from openpnm.algorithms._solution import TransientSolution
from openpnm.integrators import Integrator

__all__ = ['ScipyRK45']


class ScipyRK45(Integrator):
    """Integrator class based on SciPy's implementation of RK45"""

    def __init__(self, atol=1e-6, rtol=1e-6, verbose=False, linsolver=None):
        self.atol = atol
        self.rtol = rtol
        self.verbose = verbose
        self.linsolver = linsolver

    def solve(self, rhs, x0, tspan, saveat, **kwargs):
        """
        Solves the system of ODEs defined by dy/dt = rhs(t, y).

        Parameters
        ----------
        rhs : function handle
            RHS vector in the system of ODEs defined by dy/dt = rhs(t, y)
        x0 : array_like
            Initial value for the system of ODEs
        tspan : array_like
            2-element tuple (or array) representing the timespan for the
            system of ODEs
        saveat : float or array_like
            If float, defines the time interval at which the solution is
            to be stored. If array_like, defines the time points at which
            the solution is to be stored.
        **kwargs : keyword arguments
            Other keyword arguments that might get used by the integrator

        Returns
        -------
        TransientSolution
            Solution of the system of ODEs stored in a subclass of numpy's
            ndarray with some added functionalities (ex. you can get the
            solution at intermediate time points via: y = soln(t_i)).

        """
        options = {
            "atol": self.atol,
            "rtol": self.rtol,
            "t_eval": saveat,
            # FIXME: uncomment next line when/if scipy#11815 is merged
            # "verbose": self.verbose,
        }
        sol = solve_ivp(rhs, tspan, x0, method="RK45", **options)
        if sol.success:
            return TransientSolution(sol.t, sol.y)
        raise Exception(sol.message)
