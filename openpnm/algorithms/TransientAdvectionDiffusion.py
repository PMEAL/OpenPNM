from openpnm.algorithms import TransientReactiveTransport, AdvectionDiffusion
from openpnm.core import logging
logger = logging.getLogger(__name__)


class TransientAdvectionDiffusion(TransientReactiveTransport,
                                  AdvectionDiffusion):
    r"""
    A subclass of GenericTransport to perform steady and transient simulations
    of pure diffusion and advection diffusion problems.

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        # Apply any received settings to overwrite defaults
        self.settings.update(settings)
        # Save A matrix of the steady sys of eqs
        self.A = self._build_A()
        # Define _coef as the molar density
        phase = self.project.phases()[self.settings['phase']]
        self._coef = phase[self.settings['molar_density']]
