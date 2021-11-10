from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging, Docorator, GenericSettings
logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.get_sections(base='FourierConductionSettings',
                     sections=['Parameters'])
@docstr.dedent
class FourierConductionSettings(GenericSettings):
    r"""

    Parameters
    ----------
    %(GenericTransportSettings.parameters)s
    quantity : str (default = ``'pore.temperature'``
        The name of the physical quantity to be calculated
    conductance : str (default = ``'pore.thermal_conductance'``)
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
    quantity = 'pore.temperature'
    conductance = 'throat.thermal_conductance'


class FourierConduction(ReactiveTransport):
    r"""
    A subclass of GenericLinearTransport to simulate heat conduction.

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update(settings)
        self.settings._update_settings_and_docs(FourierConductionSettings())
