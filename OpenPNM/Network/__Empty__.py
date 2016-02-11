"""
===============================================================================
Empty: Generate an empty Network with a predefined number of pores and throats
===============================================================================

"""

import scipy as sp
from OpenPNM.Network import GenericNetwork
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class Empty(GenericNetwork):
    r"""
    Creates an empty Network but with a predefined number of pores and throats.
    This allows all of the data assignment to function properly, meaning that
    additional data can be added using the standard dictionary syntax
    (e.g. ``net['pore.data']``).  This class is intended for custom built
    networks where data is added incrementally.

    """

    def __init__(self, Np, Nt, **kwargs):
        super().__init__(**kwargs)
        self.update({'pore.all': sp.ones((Np,), dtype=bool)})
        self.update({'throat.all': sp.ones((Nt,), dtype=bool)})
