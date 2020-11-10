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
            self.setup(phase=phase)

    def setup(self, phase=None, quantity='', conductance='',
              diffusive_conductance='', ion='', t_initial=None, t_final=None,
              t_step=None, t_output=None, t_tolerance=None, t_precision=None,
              t_scheme='', **kwargs):
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        if diffusive_conductance:
            self.settings['diffusive_conductance'] = diffusive_conductance
        if ion:
            self.settings['ion'] = ion
        if t_initial is not None:
            self.settings['t_initial'] = t_initial
        if t_final is not None:
            self.settings['t_final'] = t_final
        if t_step is not None:
            self.settings['t_step'] = t_step
        if t_output is not None:
            self.settings['t_output'] = t_output
        if t_tolerance is not None:
            self.settings['t_tolerance'] = t_tolerance
        if t_precision is not None:
            self.settings['t_precision'] = t_precision
        if t_scheme:
            self.settings['t_scheme'] = t_scheme
        self.settings.update(kwargs)
