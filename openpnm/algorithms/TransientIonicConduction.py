from openpnm.algorithms import TransientReactiveTransport, IonicConduction
from openpnm.utils import logging, Docorator, GenericSettings
logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.get_sections(base='TransientIonicConductionSettings',
                     sections=['Parameters'])
@docstr.dedent
class TransientIonicConductionSettings(GenericSettings):
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
        super().__init__(**kwargs)
        c = TransientIonicConductionSettings()
        self.settings._update_settings_and_docs(c)
        self.settings.update(settings)
        if phase is not None:
            self.settings['phase'] = phase.name
