from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging, Docorator, SettingsAttr
logger = logging.getLogger(__name__)
docstr = Docorator()


class FourierConductionSettings:
    quantity = 'pore.temperature'
    conductance = 'throat.thermal_conductance'


class FourierConduction(ReactiveTransport):
    r"""
    A subclass of GenericLinearTransport to simulate heat conduction.

    """

    def __init__(self, settings={}, **kwargs):
        self.settings = SettingsAttr(FourierConductionSettings, settings)
        super().__init__(settings=self.settings, **kwargs)
