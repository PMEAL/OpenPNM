"""
module __TestNet__: Generate simple cubic networks
==========================================================

.. warning:: The classes of this module should be loaded through the 'Topology.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM

import scipy as sp
from OpenPNM.Network import GenericNetwork

class TestNet(GenericNetwork):
    r"""
    A small cubic network for quick testing purposes

    Parameters
    ----------
    This class accepts no arguments

    """

    def __init__(self, name='test_net', **kwargs):
        super(TestNet, self).__init__(name, **kwargs)
        self.generate()

    def generate(self):
        '''
        Create test network, of cubic geometry [5,5,5]

        Parameters
        ----------
        This network type accepts no arguments
        '''
        self._generate_setup()
        self._generate_pores()
        self._generate_throats()
        self._add_labels()
        return self

    def _generate_setup(self):
        r"""
        Perform applicable preliminary checks and calculations required for generation
        """
        self._Nx = 5
        self._Ny = 5
        self._Nz = 5
        self._Lc = 1
        self._Lx = sp.float16(self._Nx*self._Lc)
        self._Ly = sp.float16(self._Ny*self._Lc)
        self._Lz = sp.float16(self._Nz*self._Lc)

    def _generate_pores(self):
        r"""
        Generate the pores (coordinates, numbering and types)
        """
        Nx = self._Nx
        Ny = self._Ny
        Nz = self._Nz
        Lc = self._Lc
        Np = Nx*Ny*Nz
        ind = sp.arange(0,Np)
        self.set_pore_data(prop='numbering',data=ind)
        self.set_pore_info(label='all',locations=sp.ones_like(ind))
        pore_coords = Lc/2+Lc*sp.array(sp.unravel_index(ind, dims=(Nx, Ny, Nz), order='F'),dtype=sp.float64).T
        self.set_pore_data(prop='coords',data=pore_coords)

    def _generate_throats(self):
        r"""
        Generate the throats (connections, numbering and types)
        """
        Nx = self._Nx
        Ny = self._Ny
        Nz = self._Nz
        Np = Nx*Ny*Nz
        ind = sp.arange(0,Np)
        #Generate throats based on pattern of the adjacency matrix
        tpore1_1 = ind[(ind%Nx)<(Nx-1)]
        tpore2_1 = tpore1_1 + 1
        tpore1_2 = ind[(ind%(Nx*Ny))<(Nx*(Ny-1))]
        tpore2_2 = tpore1_2 + Nx
        tpore1_3 = ind[(ind%Np)<(Nx*Ny*(Nz-1))]
        tpore2_3 = tpore1_3 + Nx*Ny
        tpore1 = sp.hstack((tpore1_1,tpore1_2,tpore1_3))
        tpore2 = sp.hstack((tpore2_1,tpore2_2,tpore2_3))
        connections = sp.vstack((tpore1,tpore2)).T
        connections = connections[sp.lexsort((connections[:, 1], connections[:, 0]))]
        self.set_throat_data(prop='numbering',data=sp.arange(0,sp.shape(tpore1)[0]))
        self.set_throat_info(label='all',locations=sp.ones_like(sp.arange(0,sp.shape(tpore1)[0])))
        self.set_throat_data(prop='conns',data=connections)       
        
    def _add_labels(self):
        coords = self.get_pore_data(prop='coords')
        self.set_pore_info(label='front',locations=coords[:,0]<=self._Lc)
        self.set_pore_info(label='left',locations=coords[:,1]<=self._Lc)
        self.set_pore_info(label='bottom',locations=coords[:,2]<=self._Lc)
        self.set_pore_info(label='back',locations=coords[:,0]>=(self._Lc*(self._Nx-1)))
        self.set_pore_info(label='right',locations=coords[:,1]>=(self._Lc*(self._Ny-1)))
        self.set_pore_info(label='top',locations=coords[:,2]>=(self._Lc*(self._Nz-1)))
        self.set_pore_info(label='internal',locations=self.get_pore_indices())
        for item in ['top','bottom','left','right','front','back']:
            ps = self.get_pore_indices(item)
            ts = self.find_neighbor_throats(ps)
            ps = self.find_connected_pores(ts)
            ps0 = self.get_pore_info(label=item)[ps[:,0]]
            ps1 = self.get_pore_info(label=item)[ps[:,1]]
            ts = ts[ps1*ps0]
            self.set_throat_info(label=item,locations=ts)
        self.set_throat_info(label='internal',locations=self.get_throat_indices())

if __name__ == '__main__':
    pn = OpenPNM.Network.TestNet()
    print(pn.name)
