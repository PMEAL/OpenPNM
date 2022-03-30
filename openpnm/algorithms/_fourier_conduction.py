import logging
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import Docorator, SettingsAttr
logger = logging.getLogger(__name__)
docstr = Docorator()

__all__ = ['FourierConduction']


class FourierConductionSettings:
    quantity = 'pore.temperature'
    conductance = 'throat.thermal_conductance'


class FourierConduction(ReactiveTransport):
    r"""
    A subclass of GenericLinearTransport to simulate heat conduction.

    """

    def __init__(self, settings=None, **kwargs):
        self.settings = SettingsAttr(FourierConductionSettings, settings)
        super().__init__(settings=self.settings, **kwargs)
