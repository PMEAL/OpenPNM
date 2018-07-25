from openpnm.algorithms import TransientReactiveTransport, FickianDiffusion
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class TransientFickianDiffusion(TransientReactiveTransport, FickianDiffusion):
    r"""
    A subclass of GenericTransport to simulate diffusion.

    """

    def __init__(self, settings={}, **kwargs):
        def_set = {'gui': {'setup':        {'quantity': '',
                                            'conductance': '',
                                            't_initial': None,
                                            't_final': None,
                                            't_step': None,
                                            't_output': None,
                                            't_tolerance': None,
                                            't_scheme': ''
                                            },
                           'set_IC':       {'values': None,
                                            },
                           'set_rate_BC':  {'pores': None,
                                            'values': None,
                                            },
                           'set_value_BC': {'pores': None,
                                            'values': None},
                           'set_source':   {'pores': None,
                                            'propname': '',
                                            },
                           }
                   }
        super().__init__(**kwargs)
        self.settings.update(def_set)
        self.settings.update(settings)

    def setup(self, phase=None, quantity='', conductance='',
              t_initial=None, t_final=None, t_step=None, t_output=None,
              t_tolerance=None, t_scheme='', **kwargs):
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        if t_initial:
            self.settings['t_initial'] = t_initial
        if t_final:
            self.settings['t_final'] = t_final
        if t_step:
            self.settings['t_step'] = t_step
        if t_output:
            self.settings['t_output'] = t_output
        if t_tolerance:
            self.settings['t_tolerance'] = t_tolerance
        if t_scheme:
            self.settings['t_scheme'] = t_scheme
        self.settings.update(kwargs)
