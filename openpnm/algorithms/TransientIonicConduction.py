from copy import deepcopy
from openpnm.algorithms import TransientReactiveTransport, IonicConduction
from openpnm.utils import logging, Docorator, SettingsAttr
logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.dedent
class TransientIonicConductionSettings:
    r"""

    Parameters
    ----------
    ##

    Other Parameters
    ----------------

    **The following parameters pertain to steady-state IonicConduction**

    %(IonicConductionSettings.parameters)s

    ----

    **The following parameters pertain to the ReactiveTransport class**

    %(ReactiveTransportSettings.other_parameters)s

    ----

    **The following parameters pertain to the GenericTransport class**

    %(GenericTransportSettings.other_parameters)s

    """
    quantity = 'pore.potential'
    conductance = 'throat.ionic_conductance'
    charge_conservation = 'electroneutrality'
    cache_A = False
    cache_b = False


class TransientIonicConduction(TransientReactiveTransport,
                               IonicConduction):
    r"""
    A subclass of GenericTransport to perform steady and transient simulations
    of pure diffusion and advection-diffusion problems.

    """

    def __init__(self, settings={}, phase=None, **kwargs):
        self.settings = SettingsAttr(TransientIonicConductionSettings, settings)
        super().__init__(settings=self.settings, **kwargs)
        if phase is not None:
            self.settings['phase'] = phase.name
