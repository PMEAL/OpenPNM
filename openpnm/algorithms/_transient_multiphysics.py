import numpy as np
from openpnm.utils import logging, SettingsAttr
from openpnm.integrators import ScipyRK45
from openpnm.algorithms import TransientReactiveTransport
from openpnm.algorithms import GenericAlgorithm
logger = logging.getLogger(__name__)

class TransientMultiphysicsSettings:
    r"""

    Parameters
    ----------
    ...

    """
    algorithms = []

class TransientMultiphysics(GenericAlgorithm):
    r"""
    
    A class for transient multiphysics simulations
    
    """
    
    def __init__(self, algorithms, settings=None, **kwargs):
        self.algorithms = algorithms # temporarily until I get settings to work
        self.settings = SettingsAttr(TransientMultiphysicsSettings, settings)
        
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
            solution object).

        """
        logger.info('Running TransientMultiphysics')
        if np.isscalar(saveat):
            saveat = np.arange(*tspan, saveat)
        if (saveat is not None) and (tspan[1] not in saveat):
            saveat = np.hstack((saveat, [tspan[1]]))
        integrator = ScipyRK45() if integrator is None else integrator
        # for each algorithm
        algs = self.algorithms
        for i, alg in enumerate(algs):            
            # Perform pre-solve validations
            alg._validate_settings()
            alg._validate_data_health()
            # Write x0 to algorithm the obj (needed by _update_iterative_props)
            x0_i = self._get_x0(x0, i) # write helper method _get_x0
            alg['pore.ic'] = x0_i = np.ones(alg.Np, dtype=float) * x0_i
            alg._merge_inital_and_boundary_values()
        # Build RHS (dx/dt = RHS), then integrate the system of ODEs
        rhs = self._build_rhs(algs)
        # Integrate RHS using the given solver
        self.soln = integrator.solve(rhs, x0, tspan, saveat)
        # add solution to each algorithm
        for i, alg in enumerate(algs):
            alg.soln = self.soln[i*alg.Np:(i+1)*alg.Np,:]
        return self.soln

    def _run_special(self, x0): ...

    def _build_rhs(self, algs):

        def ode_func(t, y):
            # initialize rhs
            rhs = []
            for i, alg in enumerate(algs):
                # get x from y, assume alg.Np is same for all algs
                x = self._get_x0(y, i) # again use helper function
                # store x onto algorithm,
                alg.x = x
                # alg[alg.settings["quantity"]] = x
                # build A and b
                alg._update_A_and_b()
                A = alg.A.tocsc()
                b = alg.b
                # retrieve volume
                V = alg.network[alg.settings["pore_volume"]]
                # calcualte rhs
                rhs_alg = np.hstack(-A.dot(x) + b)/V
                rhs = np.hstack((rhs, rhs_alg))    
            return rhs

        return ode_func
    
    def _get_x0(self, x0, i):
        Np = self.algorithms[0].Np # fix when we change settings
        x0 = x0[i*Np:(i+1)*Np] 
        return x0