"""
module __GenericVisualization__: Base class to visualize networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Visualization/__init__.py' file.

"""

import OpenPNM
import scipy as sp

class GenericVisualization(OpenPNM.Utilities.OpenPNMbase):
    r"""
    GenericVisualization - Base class to visualize networks    
    
    Parameters
    ----------
    
    Examples
    --------
    >>> print 'nothing yet'
    
    .. note:: 
    n/a
    
    """
    
    def __init__(self,net=OpenPNM.Network.GenericNetwork,**kwords):
        r"""
        Initialize
        """
        super(GenericVisualization,self).__init__(**kwords)
        self._logger.debug("Construct class")
        self._net = net
        
if __name__ =="__main__":
    test = GenericVisualization(loggername="TestGenericVis")