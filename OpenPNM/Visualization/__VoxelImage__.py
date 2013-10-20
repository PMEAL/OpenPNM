"""
module __GenericVisualization__: Base class to visualize networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Visualization/__init__.py' file.

"""

import OpenPNM
import scipy as sp

from __GenericVisualization__ import GenericVisualization

class VoxelImage(GenericVisualization):
    r"""
    Generate a 3D voxel image of the network

    Parameters
    ----------

    Examples
    --------
    >>> print 'nothing yet'

    Notes
    -----



    """

    def __init__(self,**kwargs):
        r"""
        Initialize
        """
        super(VoxelImage,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")

    def write(self):
        r"""
        Not implemented yet

        Parameters
        ----------
        net : OpenPNM Network Object

        filename : string
            Full path to desired file location

        scaling_factor : int, optional
            Not sure what this does
        """

        print 'nothing yet'


