import logging
from openpnm.algorithms import TransientReactiveTransport, FickianDiffusion
logger = logging.getLogger(__name__)

__all__ = ['TransientFickianDiffusion']


class TransientFickianDiffusion(TransientReactiveTransport, FickianDiffusion):
    r"""
    A class to simulate transient diffusion with reactions
    """

    def __init__(self, **kwargs):
        if 'name' not in kwargs.keys():
            kwargs['name'] = 'trans_fick_01'
        super().__init__(**kwargs)
