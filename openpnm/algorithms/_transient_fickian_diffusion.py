import logging
from openpnm.algorithms import TransientReactiveTransport, FickianDiffusion


__all__ = ['TransientFickianDiffusion']


logger = logging.getLogger(__name__)


class TransientFickianDiffusion(TransientReactiveTransport, FickianDiffusion):
    r"""
    A class to simulate transient diffusion with reactions
    """

    def __init__(self, name='trans_fick_?', **kwargs):
        super().__init__(name=name, **kwargs)
