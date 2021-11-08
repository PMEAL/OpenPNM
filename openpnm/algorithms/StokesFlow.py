import numpy as np
import scipy as sp
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging, GenericSettings, Docorator
logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.get_sections(base='StokesFlowSettings',
                     sections=['Parameters'])
@docstr.dedent
class StokesFlowSettings(GenericSettings):
    r"""

    Parameters
    ----------
    %(GenericTransportSettings.parameters)s
    quantity : str (default = 'pore.pressure')
        The name of the physical quantity to be calculated
    conductance : str (default = 'throat.hydraulic_conductance')
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
    quantity = 'pore.pressure'
    conductance = 'throat.hydraulic_conductance'


class StokesFlow(ReactiveTransport):
    r"""
    A subclass of GenericLinearTransport to simulate viscous flow.

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings._update_settings_and_docs(StokesFlowSettings())
        self.settings.update(settings)