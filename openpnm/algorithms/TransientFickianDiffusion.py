from openpnm.algorithms import TransientReactiveTransport, FickianDiffusion
from openpnm.utils import logging, Docorator
logger = logging.getLogger(__name__)
docstr = Docorator()


class TransientFickianDiffusion(TransientReactiveTransport, FickianDiffusion):
    r"""
    A class to simulate transient diffusion with reactions

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update(settings)  # Add user supplied settings
