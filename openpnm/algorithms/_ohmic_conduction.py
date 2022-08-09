import logging
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import Docorator


docstr = Docorator()
logger = logging.getLogger(__name__)


__all__ = ['OhmicConduction']


@docstr.dedent
class OhmicConductionSettings:
    r"""

    Parameters
    ----------
    %(ReactiveTransportSettings.parameters)s

    """
    quantity = 'pore.voltage'
    conductance = 'throat.electrical_conductance'


class OhmicConduction(ReactiveTransport):
    r"""
    A subclass of LinearTransport to simulate electron and ionic
    conduction.  The 2 main roles of this subclass are to set the default
    property names and to implement a method for calculating the effective
    conductivity of the network.

    """

    def __init__(self, name='ohmic_?', **kwargs):
        super().__init__(name=name, **kwargs)
        self.settings._update(OhmicConductionSettings())
