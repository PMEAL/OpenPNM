from openpnm.algorithms import TransientReactiveTransport, NernstPlanck
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class TransientNernstPlanck(TransientReactiveTransport, NernstPlanck):
    r"""
    A subclass of GenericTransport to perform steady and transient simulations
    of pure diffusion, advection-diffusion and advection-diffusion with
    migration.

    """

    def __init__(self, settings={}, phase=None, ion='', **kwargs):
        def_set = {'phase': None,
                   'quantity': 'pore.concentration.'+ion,
                   'conductance': 'throat.ad_dif_mig_conductance.'+ion,
                   'diffusive_conductance': 'throat.diffusive_conductance.' +
                   ion,
                   'ion': ion,
                   'gui': {'setup':        {'phase': None,
                                            'quantity': '',
                                            'conductance': '',
                                            'diffusive_conductance': '',
                                            'ion': '',
                                            't_initial': None,
                                            't_final': None,
                                            't_step': None,
                                            't_output': None,
                                            't_tolerance': None,
                                            't_scheme': ''},
                           'set_IC':       {'values': None},
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
