import numpy as np
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging, GenericSettings, Docorator
from openpnm.solvers import ScipyRK45
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.get_sections(base='TransientReactiveTransportSettings',
                     sections=['Parameters', 'Other Parameters'])
@docstr.dedent
class TransientReactiveTransportSettings(GenericSettings):
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

    phase = None
    saveat = None
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
    - Either a Network or a Project must be specified.

    """

    def __init__(self, settings={}, phase=None, **kwargs):
        super().__init__(**kwargs)
        self.settings._update_settings_and_docs(TransientReactiveTransportSettings)
        self.settings.update(settings)
        if phase is not None:
            self.settings['phase'] = phase.name
        self["pore.ic"] = np.nan

    def run(self, x0, tspan, saveat=None, solver=None):
        """
        Parameters
        ----------
        x0 : ndarray or float
            Array (or scalar) containing initial condition values.

        """
        logger.info('Running TransientTransport')
        solver = ScipyRK45() if solver is None else solver
        # Perform pre-solve validations
        self._validate_settings()
        self._validate_data_health()
        # Write x0 to algorithm the obj (needed by _update_iterative_props)
        x0 = np.ones(self.Np) * x0
        self._overwrite_ICs_with_BCs()
        # Build RHS (dx/dt = RHS), then integrate the system of ODEs
        rhs = self._build_rhs()
        # Integrate RHS using the given solver
        return solver.solve(rhs, x0, tspan, saveat)

    def _run_special(self, x0): ...

    def _build_rhs(self):
        # TODO: add a cache mechanism
        self._build_A()
        self._build_b()
        self._apply_BCs()
        self._apply_sources()
        A = self.A.tocsc()
        b = self.b
        V = self.network[self.settings["pore_volume"]]

        def ode_func(t, y, A, b, V):
            return (-A.dot(y) + b) / V  # much faster than A*y

        return lambda t, y: ode_func(t, y, A, b, V)

    def _overwrite_ICs_with_BCs(self):
        ic_vals = self['pore.ic']
        # Ensure the given initial conditions have any value BC inserted
        bc_pores = ~np.isnan(self['pore.bc_value'])
        ic_vals[bc_pores] = self['pore.bc_value'][bc_pores]
        # Write values to self to to quantity, ic and t=0 array
        quantity = self.settings['quantity']
        self[quantity] = ic_vals
