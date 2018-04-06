from openpnm.algorithms import GenericTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class OhmicConduction(GenericTransport):
    r"""
    A subclass of GenericLinearTransport to simulate electron and ionic
    conduction.  The 2 main roles of this subclass are to set the default
    property names and to implement a method for calculating the effective
    conductivity of the network.

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update({'quantity': 'pore.voltage',
                              'conductance': 'throat.electrical_conductance'})
        self.settings.update(settings)

    def calc_effective_conductivity(self):
        r"""
        This calculates the effective electrical conductivity.
        """
        return self._calc_eff_prop()
