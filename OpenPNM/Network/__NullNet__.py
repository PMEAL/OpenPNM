"""
module __NullNet__: Generate empty networks
==========================================================

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM

import scipy as sp
from OpenPNM.Network import GenericNetwork


class NullNet(GenericNetwork):
    r"""
    An empty network for testing purposes

    Parameters
    ----------
    name : string, required
        A unique name to identify object

    """

    def __init__(self, **kwargs):
        super(NullNet, self).__init__(**kwargs)

    def generate(self, **params):
        '''
        Create empty network with no pores or throats for testing purposes

        Parameters
        ----------
        This network type accepts no arguments
        '''
        self._generate_setup(**params)
        self._add_pores()
        self._add_throats()
        return self

    def _generate_setup(self,
                        Nt = 10,
                        Np = 10,
                        **params):
        r"""
        """
        self._Np = params['Np']
        self._Nt = params['Nt']
        
    def _add_pores(self):
        self.add_pore_info(label='all', locations=sp.ones((self._Np,), dtype=bool))
        
    def _add_throats(self):
        self.add_throat_info(label='all', locations=sp.ones((self._Nt,), dtype=bool))

if __name__ == '__main__':
    pn = OpenPNM.Network.NullNet(name='null').generate(Np=10,Nt=10)
    print(pn.name)
