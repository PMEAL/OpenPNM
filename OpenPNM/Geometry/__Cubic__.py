"""
module __Cubic__: Generate simple cubic networks
==========================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""
#Test
import OpenPNM
import scipy as sp
import numpy as np
import scipy.stats as spst
import scipy.sparse as sprs
import scipy.spatial as sptl
import scipy.ndimage as spim
import matplotlib.pyplot as plt
import itertools as itr

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
        Lc = self._Lc
        params['divisions'] = [self._Nx, self._Ny, self._Nz]
        temp_divisions = sp.copy(params['divisions']).tolist()
        
        for i in range(0,6):
            params['divisions'] = temp_divisions
            temp_divisions = sp.copy(params['divisions']).tolist()
            params['divisions'][i%3] = 1
            edge_net = super(Cubic, self).generate(**params)
            displacement = [0,0,0]
            if i < 3:
                displacement[i%3] = net.pore_properties['coords'][:,i%3].max() + 0.5*Lc
            else:
                displacement[i%3] = -1*(net.pore_properties['coords'][:,i%3].min() + 0.5*Lc)
            self.translate_coordinates(edge_net,displacement)
            self.stitch(net,edge_net,edge = i+1)
            

        return net    
        self._logger.debug("generate_boundaries: End of method")

    def stitch(self,net1,net2,edge = 0):
        r"""
        Stitch two networks together
        
        Parameters
        ----------
        net1 : OpenPNM Network Object
            The network that is stiched to
        
        net2 : OpenPNM Network Object
            The network that is stitched
        
        """
        net2.pore_properties['numbering']   = len(net1.pore_properties['numbering']) + net2.pore_properties['numbering']
        net1.pore_properties['numbering']   = sp.concatenate((net1.pore_properties['numbering'],net2.pore_properties['numbering']),axis=0)        
        net1.pore_properties['volume']      = sp.concatenate((net1.pore_properties['volume'],net2.pore_properties['volume']),axis = 0)
        net1.pore_properties['seed']        = sp.concatenate((net1.pore_properties['seed'],net2.pore_properties['seed']),axis = 0)
        net2.pore_properties['type']        = sp.repeat(edge,len(net2.pore_properties['type']))
        net1.pore_properties['type']        = sp.concatenate((net1.pore_properties['type'],net2.pore_properties['type']),axis = 0)
        net1.pore_properties['diameter']    = sp.concatenate((net1.pore_properties['diameter'],net2.pore_properties['diameter']),axis = 0)
        net1.pore_properties['coords']      = sp.concatenate((net1.pore_properties['coords'],net2.pore_properties['coords']),axis = 0)
        
        net2.throat_properties['numbering'] = len(net1.throat_properties['numbering']) + net2.throat_properties['numbering']
        net1.throat_properties['numbering'] = sp.concatenate((net1.throat_properties['numbering'],net2.throat_properties['numbering']),axis=0)
        net1.throat_properties['volume']    = sp.concatenate((net1.throat_properties['volume'],net2.throat_properties['volume']),axis=0)
        net1.throat_properties['diameter']  = sp.concatenate((net1.throat_properties['diameter'],net2.throat_properties['diameter']),axis=0)
        net1.throat_properties['length']    = sp.concatenate((net1.throat_properties['length'],net2.throat_properties['length']),axis=0)
        net1.throat_properties['seed']      = sp.concatenate((net1.throat_properties['seed'],net2.throat_properties['seed']),axis=0)
        net2.throat_properties['type']      = sp.repeat(edge,len(net2.throat_properties['type']))
        net1.throat_properties['type']      = sp.concatenate((net1.throat_properties['type'],net2.throat_properties['type']),axis=0)
        
        pts = net1.pore_properties['coords']
        tri = sptl.Delaunay(pts)
        adj_mat = (sp.zeros((len(pts),len(pts)))-1).copy()
        dist_comb = list(itr.combinations_with_replacement(range(4),2))
        
        for i in range(len(pts)):
            for j in range(len(dist_comb)):
                point_1 = tri.simplices[i,dist_comb[j][0]]
                point_2 = tri.simplices[i,dist_comb[j][1]]
                coords_1 = tri.points[point_1]
                coords_2 = tri.points[point_2]
                adj_mat[point_1,point_2] = self._net.fastest_calc_dist(coords_1,coords_2)

        print adj_mat
        #plt.triplot(pts[:,0], pts[:,1], tri.simplices.copy())
        #plt.plot(pts[:,0], pts[:,1], 'o')
        #plt.show()


if __name__ == '__main__':
    test=Cubic(lattice_spacing=1.0,domain_size=[3,3,3],loggername='TestCubic')
    pn = test.generate()