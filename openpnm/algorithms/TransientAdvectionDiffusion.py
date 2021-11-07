from openpnm.algorithms import TransientReactiveTransport, AdvectionDiffusion
from openpnm.utils import logging, Docorator
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.dedent
class TransientAdvectionDiffusionSettings:
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
        self.settings._update(TransientAdvectionDiffusionSettings, docs=True)
        self.settings._update(settings)
        if phase is not None:
            self.settings['phase'] = phase.name
