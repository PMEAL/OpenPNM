import numpy as np
from copy import deepcopy
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging, Docorator
from openpnm.integrators import ScipyRK45
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

    quantity : str
        The name of the physical quantity to be calculated
    conductance : str
        The name of the pore-scale transport conductance values. These are
        typically calculated by a model attached to a *Physics* object
        associated with the given *Phase*.
    pore_volume : str
        The name of the pore volume property to use in setting up the transient
        system. Default is 'pore.volume' but 'pore.volume_effective' could be
        used if needed.

    Other Parameters
    ----------------
    saveat : float, array_like
        List of time points at which the results are to be stored (if
        array_like is passed), OR the time interval (if float is passed).

    ----

    **The following parameters pertain to the ReactiveTransport class**

    %(ReactiveTransportSettings.other_parameters)s

    ----

    **The following parameters pertain to the GenericTransport class**

    %(GenericTransportSettings.other_parameters)s

    """
    pore_volume = 'pore.volume'


class TransientReactiveTransport(ReactiveTransport):
    r"""
    A subclass of ReactiveTransport for transient simulations.

    Parameters
    ----------
    network : GenericNetwork
        The Network with which this algorithm is associated.
    project : Project
        The Project with which this algorithm is associated.

    Notes
    -----
    Either a Network or a Project must be specified.

    """

    def __init__(self, settings={}, phase=None, **kwargs):
        self.settings._update(TransientReactiveTransportSettings, docs=True)
        self.settings._update(settings)  # Add user supplied settings
        super().__init__(settings=deepcopy(self.settings), **kwargs)
        if phase is not None:
            self.settings['phase'] = phase.name
        self["pore.ic"] = np.nan

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
        if (saveat is not None) and (tspan[1] not in saveat):
            saveat = np.hstack((saveat, [tspan[1]]))
        integrator = ScipyRK45() if integrator is None else integrator
        # Perform pre-solve validations
        self._validate_settings()
        self._validate_data_health()
        # Write x0 to algorithm the obj (needed by _update_iterative_props)
        self['pore.ic'] = x0 = np.ones(self.Np, dtype=float) * x0
        self._merge_inital_and_boundary_values()
        # Build RHS (dx/dt = RHS), then integrate the system of ODEs
        rhs = self._build_rhs()
        # Integrate RHS using the given solver
        self.soln = integrator.solve(rhs, x0, tspan, saveat)
        return self.soln

    def _run_special(self, x0): ...

    def _build_rhs(self):

        def ode_func(t, y):
            # TODO: add a cache mechanism
            self[self.settings["quantity"]] = y
            TransientReactiveTransport._update_A_and_b(self)
            A = self.A.tocsc()
            b = self.b
            V = self.network[self.settings["pore_volume"]]
            return (-A.dot(y) + b) / V  # much faster than A*y

        return ode_func

    def _merge_inital_and_boundary_values(self):
        x0 = self['pore.ic']
        bc_pores = ~np.isnan(self['pore.bc_value'])
        x0[bc_pores] = self['pore.bc_value'][bc_pores]
        quantity = self.settings['quantity']
        self[quantity] = x0
