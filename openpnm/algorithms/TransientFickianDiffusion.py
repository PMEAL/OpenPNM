from openpnm.algorithms import TransientReactiveTransport, FickianDiffusion
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class TransientFickianDiffusion(TransientReactiveTransport, FickianDiffusion):
    r"""
    A subclass of GenericTransport to simulate diffusion.

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update({'molar_density': 'pore.molar_density'})
        # Apply any received settings to overwrite defaults
        self.settings.update(settings)
        # Define _coef as the molar density
        phase = self.project.phases()[self.settings['phase']]
        self._coef = phase[self.settings['molar_density']]
