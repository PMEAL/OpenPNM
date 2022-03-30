import numpy as np
from scipy.interpolate import interp1d


class Solution(np.ndarray):
    """Brief description of 'Solution'"""
    pass


class SolutionContainer(dict):
    """Brief description of 'SolutionContainer'"""
    pass


class SteadyStateSolution(np.ndarray):
    r"""Brief description of 'SteadyStateSolution'"""

    def __new__(cls, x):
        obj = np.asarray(x).view(cls)
        return obj


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
