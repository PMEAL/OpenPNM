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

    def __init__(self, network=None, project=None, settings={}, **kwargs):
        self.settings.setdefault('prefix', 'alg')
        self.settings.update(settings)
        if network is not None:
            project = network.project
        super().__init__(project=project, **kwargs)

    def results(self):
        r"""
        """
        raise NotImplementedError("This method must be subclassed")
