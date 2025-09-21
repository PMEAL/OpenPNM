import logging

import numpy as np

from openpnm.algorithms import Algorithm
from openpnm.algorithms._solution import SolutionContainer, TransientSolution
from openpnm.integrators import ScipyRK45
from openpnm.utils import Docorator

logger = logging.getLogger(__name__)
docstr = Docorator()


__all__ = [
    'TransientMultiPhysics',
]


@docstr.dedent
class TransientMultiPhysicsSettings:
    r"""

    Parameters
    ----------
    %(AlgorithmSettings.parameters)s
    algorithms: list
        List of transient algorithm objects to be solved in a coupled manner

    """
    algorithms = []


@docstr.dedent
class TransientMultiPhysics(Algorithm):
    r"""
    A subclass for transient multiphysics simulations.
    """

    def __init__(self, algorithms, **kwargs):
        super().__init__(**kwargs)
        self.settings._update(TransientMultiPhysicsSettings())
        self.settings.algorithms = [alg.name for alg in algorithms]
        self._algs = algorithms

    def run(self, x0, tspan, saveat=None, integrator=None):
        """
        Runs all of the transient algorithms simultaneoulsy and returns the
        solution.

        Parameters
        ----------
        x0 : ndarray or float
            Array (or scalar) containing initial condition values.
        tspan : array_like
            Tuple (or array) containing the integration time span.
        saveat : array_like or float, optional
            If an array is passed, it signifies the time points at which
            the solution is to be stored, and if a scalar is passed, it
            refers to the interval at which the solution is to be stored.
        integrator : Integrator, optional
            Integrator object which will be used to to the time stepping.
            Can be instantiated using openpnm.integrators module.

        Returns
        -------
        TransientSolution
            The solution object, which is basically a numpy array with
            the added functionality that it can be called to return the
            solution at intermediate times (i.e., those not stored in the
            solution object). In the case of multiphysics, the solution object
            is a combined array of solutions for each physics. The solution
            for each physics is available on each algorithm object
            independently.

        """
        logger.info('Running TransientMultiphysics')
        if np.isscalar(saveat):
            saveat = np.arange(*tspan, saveat)
        if (saveat is not None) and (tspan[1] not in saveat):
            saveat = np.hstack((saveat, [tspan[1]]))
        integrator = ScipyRK45() if integrator is None else integrator
        for i, alg in enumerate(self._algs):
            # Perform pre-solve validations
            alg._validate_settings()
            alg._validate_topology_health()
            alg._validate_linear_system()
            # Write x0 to algorithm the obj (needed by _update_iterative_props)
            x0_i = self._get_x0(x0, i)
            alg['pore.ic'] = x0_i = np.ones(alg.Np, dtype=float) * x0_i
            alg._merge_inital_and_boundary_values()
        # Build RHS (dx/dt = RHS), then integrate the system of ODEs
        rhs = self._build_rhs()
        # Integrate RHS using the given solver
        soln = integrator.solve(rhs, x0, tspan, saveat)
        # Return dictionary containing solution
        self.soln = SolutionContainer()
        for i, alg in enumerate(self._algs):
            # Slice soln and attach as TransientSolution object to each alg
            t = soln.t
            x = soln[i*alg.Np:(i+1)*alg.Np, :]
            alg.soln = TransientSolution(t, x)
            # Add solution of each alg to solution dictionary
            self.soln[alg.settings['quantity']] = alg.soln

    def _run_special(self, x0): ...

    def _build_rhs(self):
        """
        Returns a function handle, which calculates dy/dt = rhs(y, t).

        Notes
        -----
        ``y`` is a composite array that contains ALL the variables that
        the multiphysics algorithm solves for, e.g., if the constituent
        algorithms are ``TransientFickianDiffusion``, and
        ``TransientFourierConduction``, ``y[0:Np-1]`` refers to the
        concentration, and ``[Np:2*Np-1]`` refers to the temperature
        values.

        """
        def ode_func(t, y):
            # Initialize RHS
            rhs = []
            for i, alg in enumerate(self._algs):
                # Get x from y, assume alg.Np is same for all algs
                x = self._get_x0(y, i)  # again use helper function
                # Store x onto algorithm,
                alg.x = x
                # Build A and b
                alg._update_A_and_b()
                A = alg.A.tocsc()
                b = alg.b
                # Retrieve volume
                V = alg.network[alg.settings["pore_volume"]]
                # Calcualte RHS
                rhs_alg = np.hstack(-A.dot(x) + b)/V
                rhs = np.hstack((rhs, rhs_alg))
            return rhs

        return ode_func

    def _get_x0(self, x0, i):
        tmp = [alg.Np for alg in self._algs]
        idx_end = np.cumsum(tmp)
        idx_start = np.hstack((0, idx_end[:-1]))
        x0 = x0[idx_start[i]:idx_end[i]]
        return x0
