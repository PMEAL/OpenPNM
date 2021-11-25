from openpnm.algorithms import TransientReactiveTransport, AdvectionDiffusion
from openpnm.utils import logging, Docorator, SettingsAttr
docstr = Docorator()
logger = logging.getLogger(__name__)

__all__ = ['TransientAdvectionDiffusion']


@docstr.dedent
class TransientAdvectionDiffusionSettings:
    r"""
    Parameters
    ----------
    %(AdvectionDiffusionSettings.parameters)s
    %(TransientReactiveTransportSettings.parameters)s

    """
    prefix = 'trans_ad'


class TransientAdvectionDiffusion(TransientReactiveTransport,
                                  AdvectionDiffusion):
    r"""
    A subclass of GenericTransport to perform steady and transient simulations
    of pure diffusion and advection-diffusion problems.
    """

    def __init__(self, phase, settings={}, **kwargs):
        self.settings = SettingsAttr(TransientAdvectionDiffusionSettings, settings)
        super().__init__(phase=phase, settings=self.settings, **kwargs)
