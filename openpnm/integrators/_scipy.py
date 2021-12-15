import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from openpnm.integrators import Integrator

__all__ = ['ScipyRK45', 'Solution', 'TransientSolution']


class ScipyRK45(Integrator):
    """Brief description of 'ScipyRK45'"""

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


class Solution(np.ndarray):
    """Brief description of 'Solution'"""
    ...


class TransientSolution(Solution):
    """Brief description of 'TransientSolution'"""

    def __new__(cls, t, x):
        obj = np.asarray(x).view(cls)
        obj.t = np.asarray(t).view(cls)
        return obj

    def _create_interpolant(self):
        self._interpolant = interp1d(self.t, self, bounds_error=True)

    def interpolate(self, t):
        """
        Interpolates solution at time 't'.

        Parameters
        ----------
        t : float
            Time at which the solution is to be interpolated

        Returns
        -------
        ndarray
            Transient solution interpolated at the given time 't'

        Notes
        -----
        't' must reside inside the 'tspan' used during integration.

        """
        # Cache interpolant to avoid overhead
        if not hasattr(self, "_interpolant"):
            self._create_interpolant()
        return self._interpolant(t)

    __call__ = interpolate
