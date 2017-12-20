"""
===============================================================================
Bravais:
===============================================================================

"""
from openpnm.network import GenericNetwork
from openpnm import topotools
from openpnm.core import logging
logger = logging.getLogger(__name__)


class Bravais(GenericNetwork):
    r"""

    """
    def __init__(self, shape, a=1, b=1, c=1, alpha=90, beta=90, gamma=90, **kwargs):

        super().__init__()
