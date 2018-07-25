from openpnm.algorithms import TransientReactiveTransport, AdvectionDiffusion
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class TransientAdvectionDiffusion(TransientReactiveTransport,
                                  AdvectionDiffusion):
    r"""
    A subclass of GenericTransport to perform steady and transient simulations
    of pure diffusion and advection diffusion problems.

    """

    def __init__(self, settings={}, **kwargs):
        def_set = {'gui': {'setup':        {'quantity': '',
                                            'conductance': '',
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
