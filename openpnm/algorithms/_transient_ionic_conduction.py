from copy import deepcopy
from openpnm.algorithms import TransientReactiveTransport, IonicConduction
from openpnm.utils import logging, Docorator, SettingsAttr
logger = logging.getLogger(__name__)
docstr = Docorator()

__all__ = ['TransientIonicConduction']


class TransientIonicConductionSettings:
    quantity = 'pore.potential'
    conductance = 'throat.ionic_conductance'
    charge_conservation = 'electroneutrality'
    cache_A = False
    cache_b = False


class TransientIonicConduction(TransientReactiveTransport,
                               IonicConduction):
    r"""
    A subclass of GenericTransport to perform steady and transient simulations
    of pure diffusion and advection-diffusion problems.
    """

    def __init__(self, phase, settings={}, **kwargs):
        self.settings = SettingsAttr(TransientIonicConductionSettings, settings)
        super().__init__(settings=self.settings, **kwargs)
        self.settings['phase'] = phase.name
