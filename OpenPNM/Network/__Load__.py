"""
module __Cubic__: Generate simple cubic networks
==========================================================

.. warning:: The classes of this module should be loaded through the 'Topology.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM

import scipy as sp
import numpy as np
import scipy.stats as spst
import scipy.spatial as sptl
import itertools as itr
from OpenPNM.Network import GenericNetwork


class Load(GenericNetwork):
    r"""
    This class contains the methods for creating a *Cubic* network topology.  
    To invoke the actual generation it is necessary to run the `generate` method.

    Parameters
    ----------
    name : string
        A unique name for the network

    loglevel : int
        Level of the logger (10=Debug, 20=Info, 30=Warning, 40=Error, 50=Critical)
        
    loggername : string
        Overwrite the name of the logger, which defaults to the class name

    Examples
    --------
    >>> pn = OpenPNM.Network.Cubic()
    >>> pn.generate(lattice_spacing=[1],divisions=[5,5,5],add_boundaries=False)

    """

    def __init__(self, filename,**kwargs):
        super(Load, self).__init__(**kwargs)
        self._logger.debug(self.__class__.__name__+": Execute constructor")
        
        temp = filename.split('.')[0]
        temp = temp+'.npz'
        a = sp.load(temp)
        for item in a:
            self.update({item:a[item]})





        
if __name__ == '__main__':
    pn = OpenPNM.Network.Cubic(name='cubic_1',loglevel=10).generate(lattice_spacing=[1.0],domain_size=[3,3,3])
    print(pn.name)
