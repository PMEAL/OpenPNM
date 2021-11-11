from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging, Docorator, SettingsAttr
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.get_sections(base='OhmicConductionSettings',
                     sections=['Parameters'])
@docstr.dedent
class OhmicConductionSettings:
    r"""

    Parameters
    ----------
    %(GenericTransportSettings.parameters)s
    quantity : str (default = ``'pore.voltage'``)
        The name of the physical quantity to be calculated
    conductance : str (default = ``'throat.electrical_conductance'``)
        The name of the pore-scale transport conductance values. These are
        typically calculated by a model attached to a *Physics* object
        associated with the given *Phase*.

    Other Parameters
    ----------------

    **The following parameters pertain to the ReactiveTransport class**

    %(ReactiveTransportSettings.other_parameters)s

    ----

    **The following parameters pertain to the GenericTransport class**

    %(GenericTransportSettings.other_parameters)s

    """
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
