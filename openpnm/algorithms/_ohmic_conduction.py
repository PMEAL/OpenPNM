from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging, Docorator, SettingsAttr
docstr = Docorator()
logger = logging.getLogger(__name__)

__all__ = ['OhmicConduction']


@docstr.get_sections(base='OhmicConductionSettings',
                     sections=['Parameters'])
@docstr.dedent
class OhmicConductionSettings:
    quantity = 'pore.voltage'
    conductance = 'throat.electrical_conductance'


class OhmicConduction(ReactiveTransport):
    r"""
    A subclass of GenericLinearTransport to simulate electron and ionic
    conduction.  The 2 main roles of this subclass are to set the default
    property names and to implement a method for calculating the effective
    conductivity of the network.

    """

    def __init__(self, settings={}, **kwargs):
        self.settings = SettingsAttr(OhmicConductionSettings, settings)
        super().__init__(settings=self.settings, **kwargs)
