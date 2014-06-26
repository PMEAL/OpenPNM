"""
module __Load__: Load saved network data
==========================================================


"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM

import scipy as sp
from OpenPNM.Network import GenericNetwork


class Load(GenericNetwork):
    r"""
    This class loads data stored in an NPZ file from a previously saved
    network

    Parameters
    ----------
    filename : string
        The file name containing the into to be imported

    name : string
        A unique name for the network

    Examples
    --------
    >>> pn = OpenPNM.Network.Load(name='testnet',filename='test.npz')
    >>> pn.name
    'testnet'
    >>> pn.num_pores()
    10400

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
    pn = OpenPNM.Network.Load('test.npz')
    print(pn.name)
