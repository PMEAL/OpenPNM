from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class NonNewtonianStokesFlow(ReactiveTransport):
    r"""
    A subclass of GenericLinearTransport to simulate viscous flow.  The 2
    main roles of this subclass are to set the default property names and to
    implement a method for calculating the hydraulic permeability of the
    network.

    """

    def __init__(self, settings={}, phase=None, **kwargs):
        def_set = {'phase': None,
                   'quantity': 'pore.pressure',
                   'conductance': 'throat.nonNewtonian_hydraulic_conductance',
                   'gui': {'setup':        {'phase': None,
                                            'quantity': '',
                                            'conductance': ''},
                           'set_rate_BC':  {'pores': None,
                                            'values': None},
                           'set_value_BC': {'pores': None,
                                            'values': None},
                           'set_source':   {'pores': None,
                                            'propname': ''}
                           }
                   }
        super().__init__(**kwargs)
        self.settings.update(def_set)
        self.settings.update(settings)
        if phase is not None:
            self.settings['phase'] = phase.name
