from openpnm.algorithms import TransientReactiveTransport, FickianDiffusion
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class TransientFickianDiffusion(TransientReactiveTransport, FickianDiffusion):
    r"""
    A subclass of GenericTransport to simulate diffusion.

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        # Apply any received settings to overwrite defaults
        self.settings.update(settings)
