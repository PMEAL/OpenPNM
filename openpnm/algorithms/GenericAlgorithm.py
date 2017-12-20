# -*- coding: utf-8 -*-
"""
===============================================================================
module __GenericAlgorithm__: Base class to build custom algorithms
==================================================================

This generic class contains the recommended methods for subclassed algorithms.
It inherits from Core, so is Python Dict with the OpenPNM data control methods.

"""
from openpnm.core import Base
from openpnm.core import logging
logger = logging.getLogger()


class GenericAlgorithm(Base):
    r"""
    GenericAlgorithm - Base class to execute algorithms

    Parameters
    ----------
    network : OpenPNM Network Object
        The network object to which this algorithm will apply.

    name : string, optional
        Name of this algorithm

    """

    def __init__(self, network, phase, **kwargs):
        super().__init__(simulation=network.simulation, **kwargs)
        logger.name = self.name

        # Initialize label 'all' in the object's own info dictionaries
        self['pore._id'] = network['pore._id']
        self['throat._id'] = network['throat._id']
        self['pore.all'] = network['pore.all']
        self['throat.all'] = network['throat.all']
