from openpnm.algorithms import TransientReactiveTransport, AdvectionDiffusion
from openpnm.utils import logging, Docorator, GenericSettings
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.dedent
class TransientAdvectionDiffusionSettings(GenericSettings):
    r"""
    Parameters
    ----------
    %(AdvectionDiffusionSettings.parameters)s

    Other Parameters
    ----------------
    %(AdvectionDiffusionSettings.other_parameters)s

    ----

    **The following parameters pertain to the TransientReactiveTransport class**

    %(TransientReactiveTransportSettings.other_parameters)s

    ----

    **The following parameters pertain to the ReactiveTransport class**

    %(ReactiveTransportSettings.other_parameters)s

    ----

    **The following parameters pertain to the GenericTransport class**

    %(GenericTransportSettings.other_parameters)s
    """


class TransientAdvectionDiffusion(TransientReactiveTransport,
                                  AdvectionDiffusion):
    r"""
    A subclass of GenericTransport to perform steady and transient simulations
    of pure diffusion and advection-diffusion problems.

    """

    def __init__(self, settings={}, phase=None, **kwargs):
        super().__init__(**kwargs)
        self.settings._update_settings_and_docs(TransientAdvectionDiffusionSettings())
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)

    def setup(self, phase=None, quantity='', conductance='',
              diffusive_conductance='', hydraulic_conductance='', pressure='',
              s_scheme='', t_initial=None, t_final=None, t_step=None,
              t_output=None, t_tolerance=None, t_precision=None, t_scheme='',
              **kwargs):
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        if diffusive_conductance:
            self.settings['diffusive_conductance'] = diffusive_conductance
        if hydraulic_conductance:
            self.settings['hydraulic_conductance'] = hydraulic_conductance
        if pressure:
            self.settings['pressure'] = pressure
        if s_scheme:
            self.settings['s_scheme'] = s_scheme
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
