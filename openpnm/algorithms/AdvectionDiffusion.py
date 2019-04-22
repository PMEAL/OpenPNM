from openpnm.algorithms.Dispersion import Dispersion
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class AdvectionDiffusion(Dispersion):
    r"""
    A subclass of GenericTransport to simulate advection-diffusion

    """

    def __init__(self, settings={}, phase=None, **kwargs):
        def_set = {'phase': None,
                   'conductance': 'throat.ad_dif_conductance',
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
            self.setup(phase=phase)
