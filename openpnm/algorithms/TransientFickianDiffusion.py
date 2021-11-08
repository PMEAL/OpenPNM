from copy import deepcopy
from openpnm.algorithms import TransientReactiveTransport, FickianDiffusion
from openpnm.utils import logging, Docorator
logger = logging.getLogger(__name__)
docstr = Docorator()


class TransientFickianDiffusion(TransientReactiveTransport, FickianDiffusion):
    r"""
    A class to simulate transient diffusion with reactions

    """

    def __init__(self, settings={}, **kwargs):
        self.settings._update(settings)  # Add user supplied settings
        super().__init__(settings=deepcopy(self.settings), **kwargs)
