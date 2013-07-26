"""
module __CubicGenerators__: Generate simple cubic networks
==========================================================

.. warning:: The classes of this module should be loaded through the 'GEN.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import numpy as np
import scipy.sparse as sprs
import scipy.stats as spst
import scipy.ndimage as spim
from time import clock

from __GenericGenerator__ import GenericGenerator

class Custom(GenericGenerator):
    r"""
    Cubic - Class to create a basic cubic network
    
    Parameters
    ----------
    
    domain_size : list with 3 float elements 
        Shape of the cube [Lx,Ly,Lz]
    lattice_spacing : list of three floats
        Spacing between pore centers in each spatial directions
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
        
    Examples
    --------
    
    .. plot::
        
       import pylab as pl
       import OpenPNM
       gen = OpenPNM.GEN.Cubic()
       net = gen.generate()
       pl.spy(net._adjmatrix)
       pl.show()
    
    TODO:
        - Check for 3D shape
        - Check for correct divisions
    """
    
    def __init__(self,  domain_size = [],
                        divisions = [],
                        lattice_spacing = [],
                        **kwargs
                      ):

        super(Cubic,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")
        self._logger.info("Import image containing custom network shape")
        
        Np = self._Nx*self._Nx*self._Nx
        Nt = 3*Np - self._Nx*self._Ny - self._Nx*self._Nz - self._Ny*self._Nz
        
        img = np.ones((50,50,20),dtype=int)
        img[25,25,10] = 0
        img = (spim.distance_transform_edt(img)<20)
        plt.imshow(img[:,:,10])
        
        #Instantiate object
        self._net=OpenPNM.NET.GenericNetwork(num_pores=Np, num_throats=Nt)
    
    def generate_pores(self):
        r"""
        Generate the pores (coordinates, numbering and types)
        """
        self._logger.info("generate_pores: Create specified number of pores")
        Nx = self._Nx
        Ny = self._Ny
        Nz = self._Nz
        Lc = self._Lc
        Np = Nx*Ny*Nz
        ind = np.arange(0,Np)
        a = np.array(np.unravel_index(ind, dims=(Nx, Ny, Nz), order='F')).T
        self._net.pore_properties['coords'] = Lc*(0.5 + np.array(np.unravel_index(ind, dims=(Nx, Ny, Nz), order='F')).T)
        self._net.pore_properties['numbering'] = ind
        self._net.pore_properties['type']= np.zeros((Np,))
        
        self._logger.debug("generate_pores: End of method")
        
    def generate_throats(self):
        r"""
        Generate the throats (connections, numbering and types)
        """
        self._logger.info("generate_throats: Define connections between pores")
        
        Nx = self._Nx
        Ny = self._Ny
        Nz = self._Nz        
        Np = Nx*Ny*Nz
        ind = np.arange(0,Np)
        
        #Generate throats based on pattern of the adjacency matrix
        tpore1_1 = ind[(ind%Nx)<(Nx-1)]
        tpore2_1 = tpore1_1 + 1
        tpore1_2 = ind[(ind%(Nx*Ny))<(Nx*(Ny-1))]
        tpore2_2 = tpore1_2 + Nx
        tpore1_3 = ind[(ind%Np)<(Nx*Ny*(Nz-1))]
        tpore2_3 = tpore1_3 + Nx*Ny
        tpore1 = np.hstack((tpore1_1,tpore1_2,tpore1_3))
        tpore2 = np.hstack((tpore2_1,tpore2_2,tpore2_3))
        connections = np.vstack((tpore1,tpore2)).T
        connections = connections[np.lexsort((connections[:, 1], connections[:, 0]))]
        self._net.throat_properties['connections'] = connections
        self._net.throat_properties['type'] = np.zeros(np.shape(tpore1))
        self._net.throat_properties['numbering'] = np.arange(0,np.shape(tpore1)[0])
        self._logger.debug("generate_throats: End of method")
        
        
if __name__ == '__main__':
    test=Cubic(loggername='TestCubic')
    test.generate()