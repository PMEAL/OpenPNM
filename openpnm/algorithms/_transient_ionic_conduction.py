import logging
from openpnm.algorithms import TransientReactiveTransport, IonicConduction
from openpnm.utils import Docorator


__all__ = ['TransientIonicConduction']


logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.dedent
class TransientIonicConductionSettings:
    r'''

    Parameters
    ----------
    %(ReactiveTransportSettings.parameters)s
    %(TransientReactiveTransportSettings.parameters)s

    '''
    quantity = 'pore.potential'
    conductance = 'throat.ionic_conductance'
    charge_conservation = 'electroneutrality'
    cache = False


class TransientIonicConduction(TransientReactiveTransport,
                               IonicConduction):
    r"""
    A subclass of Transport to perform steady and transient simulations
    of pure diffusion and advection-diffusion problems.
    """

    def __init__(self, name='trans_ionic_#', **kwargs):
        super().__init__(name=name, **kwargs)
        self.settings._update(TransientIonicConductionSettings())
