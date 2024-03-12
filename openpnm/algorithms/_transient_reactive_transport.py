import logging
import numpy as np
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import Docorator
from openpnm.integrators import ScipyRK45
from openpnm.algorithms._solution import SolutionContainer


__all__ = ['TransientReactiveTransport']


docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.get_sections(base='TransientReactiveTransportSettings',
                     sections=['Parameters', 'Other Parameters'])
@docstr.dedent
class TransientReactiveTransportSettings:
    r"""

    Parameters
    ----------
    %(ReactiveTransportSettings.parameters)s

    """
    pore_volume = 'pore.volume'


class TransientReactiveTransport(ReactiveTransport):
    r"""
    A subclass of ReactiveTransport for transient simulations.

    Parameters
    ----------
    network : Network
        The Network with which this algorithm is associated.

    Notes
    -----
    Either a Network or a Project must be specified.

    """

    def __init__(self, phase, name='trans_react_?', **kwargs):
        super().__init__(phase=phase, name=name, **kwargs)
        self.settings._update(TransientReactiveTransportSettings())
        self.settings['phase'] = phase.name
        self["pore.ic"] = np.nan
        self._callbacks = []

    def run(self, x0, tspan, saveat=None, integrator=None):
        """
        Runs the transient algorithm and returns the solution.

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
        logger.info('Running TransientTransport')
        if np.isscalar(saveat):
            saveat = np.arange(*tspan, saveat)
        # FIXME: why do we forcibly add tspan[1] to saveat even if the user
        # didn't want to?
        if (saveat is not None) and (tspan[1] not in saveat):
            saveat = np.hstack((saveat, [tspan[1]]))
        integrator = ScipyRK45() if integrator is None else integrator
        # Perform pre-solve validations
        self._validate_settings()
        self._validate_topology_health()
        self._validate_linear_system()
        # Write x0 to algorithm the obj (needed by _update_iterative_props)
        self['pore.ic'] = x0 = np.ones(self.Np, dtype=float) * x0
        self._merge_inital_and_boundary_values()
        # Build RHS (dx/dt = RHS), then integrate the system of ODEs
        rhs = self._build_rhs()
        # Integrate RHS using the given solver
        soln = integrator.solve(rhs, x0, tspan, saveat)
        # Return solution as dictionary
        self.soln = SolutionContainer()
        self.soln[self.settings['quantity']] = soln

    def _run_special(self, x0):
        pass

    def _set_callback(self, func):
        """
        Stores the given function handle as a callback.

        Callback functions are called at the beginning of every time step.

        Parameters
        ----------
        func : function
            Function to be called as a callback. Must be a function of
            t and y.

        Returns
        -------
        None

        Examples
        --------
        >>> import openpnm as op
        >>> net = op.network.Demo([1, 2, 3])
        >>> air = op.phase.Air(network=net)
        >>> air.add_model_collection(op.models.collections.physics.standard)
        >>> trt = op.algorithms.TransientReactiveTransport(network=net, phase=air)
        >>> func = lambda t, y: print(t)
        >>> trt._set_callback(func)

        """
        if not callable(func):
            raise Exception("'func' must be a function. See Examples section.")
        self._callbacks.append(func)

    def _build_rhs(self):
        """
        Returns a function handle, which calculates dy/dt = rhs(y, t).

        Notes
        -----
        ``y`` is the variable that the algorithms solves for, e.g., for
        ``TransientFickianDiffusion``, it would be concentration.

        """
        def ode_func(t, y):
            # TODO: add a cache mechanism
            self._apply_callbacks(t, y)
            self.x = y
            self._update_A_and_b()
            A = self.A.tocsc()
            b = self.b
            V = self.network[self.settings["pore_volume"]]
            return (-A.dot(y) + b) / V  # much faster than A*y

        return ode_func

    def _merge_inital_and_boundary_values(self):
        x0 = self['pore.ic']
        bc_pores = ~np.isnan(self['pore.bc.value'])
        x0[bc_pores] = self['pore.bc.value'][bc_pores]
        quantity = self.settings['quantity']
        self[quantity] = x0

    def _apply_callbacks(self, t, y):
        for func in self._callbacks:
            func(t, y)
