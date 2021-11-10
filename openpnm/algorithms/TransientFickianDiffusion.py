from openpnm.algorithms import TransientReactiveTransport, FickianDiffusion
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class TransientFickianDiffusion(TransientReactiveTransport, FickianDiffusion):
    r"""
    A class to simulate transient diffusion with reactions

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
