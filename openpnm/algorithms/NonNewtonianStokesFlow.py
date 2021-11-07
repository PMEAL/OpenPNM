from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class NonNewtonianStokesFlowSettings:
    phase = ''
    quantity = 'pore.pressure'
    conductance = 'throat.nonNewtonian_hydraulic_conductance'


class NonNewtonianStokesFlow(ReactiveTransport):
    r"""
    A subclass of GenericLinearTransport to simulate viscous flow.  The 2
    main roles of this subclass are to set the default property names and to
    implement a method for calculating the hydraulic permeability of the
    network.

    """

    def __init__(self, settings={}, phase=None, **kwargs):
        super().__init__(**kwargs)
        self.settings._update(NonNewtonianStokesFlowSettings, docs=True)
        self.settings._update(settings)  # Add user supplied settings
        if phase is not None:
            self.settings['phase'] = phase.name
