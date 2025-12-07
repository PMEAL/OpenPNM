import logging
from openpnm.algorithms import TransientReactiveTransport, FourierConduction


__all__ = ['TransientFourierConduction']


logger = logging.getLogger(__name__)


class TransientFourierConductionSettings():
    r"""

    Parameters
    ----------
    %(TransientFourierConduction.parameters)s

    """
    quantity = 'pore.temperature'
    conductance = 'throat.thermal_conductance'
    pore_volume = 'pore.volume'


class TransientFourierConduction(TransientReactiveTransport, FourierConduction):
    r"""
    A class to simulate transient thermal diffusion with heat source
    """

    def __init__(self, name='trans_fourier_?', **kwargs):
        super().__init__(name=name, **kwargs)
