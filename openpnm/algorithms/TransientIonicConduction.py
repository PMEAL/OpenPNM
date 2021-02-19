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
            self.setup(phase=phase)

    def setup(self, phase=None, quantity='', charge_conservation=None,
              conductance='', t_initial=None, t_final=None, t_step=None,
              t_output=None, t_tolerance=None, t_precision=None, t_scheme='',
              **kwargs):
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        if charge_conservation:
            self.settings['charge_conservation'] = charge_conservation
        if t_initial is not None:
            self.settings['t_initial'] = t_initial
        if t_final is not None:
            self.settings['t_final'] = t_final
        if t_step is not None:
            self.settings['t_step'] = t_step
        if t_output is not None:
            self.settings['t_output'] = t_output
        if t_tolerance is not None:
            self.settings['t_tolerance'] = t_tolerance
        if t_precision is not None:
            self.settings['t_tolerance'] = t_precision
        if t_scheme:
            self.settings['t_scheme'] = t_scheme
        self.settings.update(kwargs)
