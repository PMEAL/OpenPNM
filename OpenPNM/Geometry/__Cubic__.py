"""
module __Cubic__: Generate simple cubic networks
==========================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import numpy as np
import scipy.stats as spst
import scipy.sparse as sprs
import scipy.spatial as sptl
import scipy.ndimage as spim

from __GenericGeometry__ import GenericGeometry

class Cubic(GenericGeometry):
    r"""
    Cubic - Class to create a basic cubic network
    
    Parameters
    ----------
    
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
        
    Examples
    --------    
    >>>print 'none yet'

    """
    
    def __init__(self, **kwargs):
        super(Cubic,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")
        
        #Instantiate pore network object
        self._net=OpenPNM.Network.GenericNetwork()

    def _generate_setup(self,   domain_size = [],
                                divisions = [3,3,3],
                                lattice_spacing = [1],
                                **params):
        r"""
        Perform applicable preliminary checks and calculations required for generation
        """
        self._logger.debug("generate_setup: Perform preliminary calculations")
        #Parse the given network size variables
        self._logger.info("Find network dimensions")
        
        self._psd = params['psd_info']
        self._tsd = params['tsd_info']
        
        if domain_size and lattice_spacing and not divisions:
            self._Lc = lattice_spacing[0]
            self._Lx = domain_size[0]
            self._Ly = domain_size[1]
            self._Lz = domain_size[2]
            self._Nx = int(self._Lx/self._Lc)
            self._Ny = int(self._Ly/self._Lc)
            self._Nz = int(self._Lz/self._Lc)
        elif divisions and lattice_spacing and not domain_size:
            self._Lc = lattice_spacing[0]
            self._Nx = divisions[0]
            self._Ny = divisions[1]
            self._Nz = divisions[2]
            self._Lx = self._Nx*self._Lc
            self._Ly = self._Ny*self._Lc
            self._Lz = self._Nz*self._Lc
        elif domain_size and divisions and not lattice_spacing:
            self._Lc = np.min(np.array(domain_size)/np.array(divisions))
            self._Nx = divisions[0]
            self._Ny = divisions[1]
            self._Nz = divisions[2]
            self._Lx = self._Nx*self._Lc
            self._Ly = self._Ny*self._Lc
            self._Lz = self._Nz*self._Lc
        else:
            self._logger.error("Exactly two of domain_size, divisions and lattice_spacing must be given")
            raise Exception('error')
        
    def _generate_pores(self):
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
        self._generate_coords(dom = (Nx,Ny,Nz), Lc = Lc)
        self._net.pore_properties['coords'] = self._temp_coords
        self._net.pore_properties['numbering'] = ind
        self._net.pore_properties['type']= np.zeros((Np,),dtype=np.int8)
        
        self._logger.debug("generate_pores: End of method")
        
    def _generate_coords(self, dom = (0,0,0), Lc = 0):
        r"""
        Generate the coordinates (coordinates, lattice Lc)
        """
        self._logger.info("generate_coordinates by dimensions and lattice parameters: Create specified number of pores")
        Np = dom[0]*dom[1]*dom[2]
        ind = np.arange(0,Np)
        self._temp_coords = Lc*(0.5 + np.array(np.unravel_index(ind, dims=(dom[0], dom[1], dom[2]), order='F')).T)
        self._logger.debug("generate_coordinates: End of method")
        
    def _generate_throats(self):
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
        self._net.throat_properties['type'] = np.zeros(np.shape(tpore1),dtype=np.int8)
        self._net.throat_properties['numbering'] = np.arange(0,np.shape(tpore1)[0])
        self._logger.debug("generate_throats: End of method")
        
    def _generate_boundaries(self,net,**params):
        self._logger.info("generate_boundaries: Define edge points of the pore network and stitch them on")
        self._generate_setup(**params)
        params['divisions'] = [self._Nx, self._Ny, self._Nz]
        
        temp_divisions = sp.copy(params['divisions']).tolist()
        params['divisions'][0] = 1
        edge_net_1 = super(Cubic, self).generate(**params)
        self.translate_coordinates(edge_net_1,[edge_net_1.pore_properties['coords'][:,0].max(),0,0])
        self.stitch(net,edge_net_1,edge = 1)
        
        params['divisions'] = temp_divisions
        temp_divisions = sp.copy(params['divisions']).tolist()
        params['divisions'][1] = 1
        edge_net_2 = super(Cubic, self).generate(**params)
        self.translate_coordinates(edge_net_2,[0,edge_net_2.pore_properties['coords'][:,1].max(),0])
        self.stitch(net,edge_net_2,edge = 2)
        
        params['divisions'] = temp_divisions
        temp_divisions = sp.copy(params['divisions']).tolist()
        params['divisions'][2] = 1
        edge_net_3 = super(Cubic, self).generate(**params)
        self.translate_coordinates(edge_net_3,[0,0,edge_net_3.pore_properties['coords'][:,2].max()])
        self.stitch(net,edge_net_3,edge = 3)
        
        temp_divisions = sp.copy(params['divisions']).tolist()
        params['divisions'][0] = 1
        edge_net_4 = super(Cubic, self).generate(**params)
        self.translate_coordinates(edge_net_1,[edge_net_4.pore_properties['coords'][:,0].max(),0,0])
        self.stitch(net,edge_net_4,edge = 4)
        
        params['divisions'] = temp_divisions
        temp_divisions = sp.copy(params['divisions']).tolist()
        params['divisions'][1] = 1
        edge_net_5 = super(Cubic, self).generate(**params)
        self.translate_coordinates(edge_net_5,[0,edge_net_5.pore_properties['coords'][:,1].max(),0])
        self.stitch(net,edge_net_5,edge = 5)
        
        params['divisions'] = temp_divisions
        temp_divisions = sp.copy(params['divisions']).tolist()
        params['divisions'][2] = 1
        edge_net_6 = super(Cubic, self).generate(**params)
        self.translate_coordinates(edge_net_6,[0,0,edge_net_6.pore_properties['coords'][:,2].max()])
        self.stitch(net,edge_net_6,edge = 6)
        return net    
        self._logger.debug("generate_boundaries: End of method")
        
if __name__ == '__main__':
    test=Cubic(lattice_spacing=1.0,domain_size=[3,3,3],loggername='TestCubic')
    pn = test.generate()