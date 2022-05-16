import logging
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import Docorator, SettingsAttr


logger = logging.getLogger(__name__)
docstr = Docorator()


__all__ = ['FourierConduction']


@docstr.dedent
class FourierConductionSettings:
    quantity = 'pore.temperature'
    conductance = 'throat.thermal_conductance'


class FourierConduction(ReactiveTransport):
    r"""
    A subclass of GenericLinearTransport to simulate heat conduction

    """

    def __init__(self, settings=None, **kwargs):
        self.settings = SettingsAttr(FourierConductionSettings, settings)
        if 'name' not in kwargs.keys():
            kwargs['name'] = 'fourier_01'
        super().__init__(settings=self.settings, **kwargs)
