# -*- coding: utf-8 -*-
"""
===============================================================================
module __GenericAlgorithm__: Base class to build custom algorithms
===============================================================================

"""
from openpnm.core import Base, logging
logger = logging.getLogger()
import numpy as np


class GenericAlgorithm(Base):
    r"""

    Parameters
    ----------
    network : OpenPNM Network Object
        The network object to which this algorithm will apply.

    name : string, optional
        Name of the algorithm

    """

    def __init__(self, network=None, project=None, settings={}, **kwargs):
        self.settings.setdefault('prefix', 'alg')
        self.settings.update(settings)
        if network is not None:
            project = network.project
        super().__init__(project=project, **kwargs)

        if project.network is not None:
            self['pore.all'] = np.ones((project.network.Np, ), dtype=bool)
            self['throat.all'] = np.ones((project.network.Nt, ), dtype=bool)

    def results(self):
        r"""
        """
        raise NotImplementedError("This method must be subclassed")
