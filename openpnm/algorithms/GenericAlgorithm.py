# -*- coding: utf-8 -*-
"""
===============================================================================
module __GenericAlgorithm__: Base class to build custom algorithms
===============================================================================

"""
from openpnm.core import Base, logging
logger = logging.getLogger()


class GenericAlgorithm(Base):
    r"""

    Parameters
    ----------
    network : OpenPNM Network Object
        The network object to which this algorithm will apply.

    name : string, optional
        Name of the algorithm

    """

    def __init__(self, network, settings={}, **kwargs):
        self.settings.setdefault('prefix', 'alg')
        self.settings.update(settings)
        super().__init__(simulation=network.simulation, **kwargs)
