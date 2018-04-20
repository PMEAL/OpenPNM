from openpnm.algorithms import FickianDiffusion
from openpnm.algorithms import TransientTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class TransientFickianDiffusion(TransientTransport, FickianDiffusion):
    r"""
    A subclass of GenericTransport to simulate diffusion.

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update({'molar_density': 'pore.molar_density'})
        # Apply any received settings to overwrite defaults
        self.settings.update(settings)
        # Save A matrix of the steady sys of eqs
        self.A = self._build_A()
        # Define _coef as the molar density
        phase = self.project.phases()[self.settings['phase']]
        self._coef = phase[self.settings['molar_density']]
