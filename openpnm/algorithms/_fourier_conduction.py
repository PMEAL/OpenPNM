import logging
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import Docorator


logger = logging.getLogger(__name__)
docstr = Docorator()


__all__ = ['FourierConduction']


@docstr.dedent
class FourierConductionSettings:
    r'''

    Parameters
    ----------
    %(ReactiveTransportSettings.parameters)s

    '''
    quantity = 'pore.temperature'
    conductance = 'throat.thermal_conductance'


class FourierConduction(ReactiveTransport):
    r"""
    A subclass of GenericLinearTransport to simulate heat conduction

    """

    def __init__(self, **kwargs):
        if 'name' not in kwargs.keys():
            kwargs['name'] = 'fourier_01'
        super().__init__(**kwargs)
        self.settings._update(FourierConductionSettings())
