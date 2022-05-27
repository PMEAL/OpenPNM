import logging
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import Docorator


__all__ = ['StokesFlow']


logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.dedent
class StokesFlowSettings:
    r'''

    Parameters
    ----------
    %(ReactiveTransportSettings.parameters)s
    '''
    quantity = 'pore.pressure'
    conductance = 'throat.hydraulic_conductance'


class StokesFlow(ReactiveTransport):
    r"""
    A subclass of GenericLinearTransport to simulate viscous flow.
    """

    def __init__(self, **kwargs):
        if 'name' not in kwargs.keys():
            kwargs['name'] = 'stokes_01'
        super().__init__(**kwargs)
        self.settings._update(StokesFlowSettings())
