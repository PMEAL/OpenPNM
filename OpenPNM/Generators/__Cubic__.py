"""
module __Cubic__: Generate simple cubic networks
==========================================================

.. warning:: The classes of this module should be loaded through the 'Generators.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import numpy as np
import scipy.sparse as sprs
import scipy.stats as spst
#import numexpr as ne
from time import clock

from __GenericGenerator__ import GenericGenerator

class Cubic(GenericGenerator):
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
    >>>print 'none yet'

    """
    
    def __init__(self,  domain_size = [],
                        divisions = [],
                        lattice_spacing = [],
                        **kwargs
                      ):
        super(Cubic,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")
        self._logger.info("Find network dimensions")
        #Parse the given network size variables
        if domain_size and lattice_spacing and not divisions:           
            self._Lc = lattice_spacing
            # Subtract all of the xyz terms for N and L by a value of 2 to compensate for boundary addition, re-add later.
            self._Nx = int(domain_size[0]/self._Lc)
            self._Ny = int(domain_size[1]/self._Lc)
            self._Nz = int(domain_size[2]/self._Lc)
            self._Lx = domain_size[0]
            self._Ly = domain_size[1]
            self._Lz = domain_size[2]
        elif divisions and lattice_spacing and not domain_size:
            print divisions
            self._Lc = lattice_spacing
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
            self._Lz = self._Nz*self._Lcp
        else:
            self._logger.error("Exactly two of domain_size, divisions and lattice_spacing must be given")
            raise Exception('error')

        Np = self._Nx*self._Ny*self._Nz
        Nt = 3*Np - self._Nx*self._Ny - self._Nx*self._Nz - self._Ny*self._Nz
        
        #Instantiate object(correct terminology?)
        self._net= OpenPNM.Network.GenericNetwork(num_pores=Np, num_throats=Nt),

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
        self._net.pore_properties['type']= np.zeros((Np,),dtype=np.int8)
        
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
        self._net.throat_properties['type'] = np.zeros(np.shape(tpore1),dtype=np.int8)
        self._net.throat_properties['numbering'] = np.arange(0,np.shape(tpore1)[0])
        self._logger.debug("generate_throats: End of method")
    
    def generate_boundary(self):
        r'''
        Generate the pores for the boundary layer. This function should be called a total of 3 times with different coordinates for x and y.
        '''
        
        # No equivalent of case switch statements. Need to rewrite more elegantly using a dictionary.
        
        if boundary[0] and boundary[1] and not boundary[2]:
            print boundary,"XYPlane"
            # This is the X and Y plane boundary
        elif boundary[0] and not boundary[1] and boundary[2]:
            print boundary,"XZPlane"
            # This is the X and Z plane boundary
        elif not boundary[0] and boundary[1] and boundary[2]:
            print boundary,"YZPlane"
            # This is the Y and Z plane boundary    

    '''        
    def add_boundaries(self):
        self._logger.debug("add_boundaries: Start of method")
        #Remove all items pertaining to previously defined boundaries (if any)
        Np = self._net.get_num_pores([0])
        Nt = self._net.get_num_throats([0])
#        for item in self._net.pore_properties.keys():
#            self._net.pore_properties[item] = self._net.pore_properties[item][0:Np]
#        for item in self._net.throat_properties.keys():
#            self._net.throat_properties[item] = self._net.throat_properties[item][0:Nt]
        pnum_orig = self._net.get_num_pores([0])
        self.add_opposing_boundaries(face=2,periodic=self._btype[0]) #x faces
        self.add_opposing_boundaries(face=3,periodic=self._btype[1]) #y faces
        self.add_opposing_boundaries(face=1,periodic=self._btype[2]) #z faces
        
        pnum_added = self._net.get_num_pores() - pnum_orig
        self._net.pore_properties['coords'] = np.concatenate((self._net.pore_properties['coords'],np.zeros((pnum_added,3))),axis=0)
        #Add 'coords' to boundaries
        #   Generate an Nx2 array, named "boundary_pore_list" that names all 
        #   pairs of pores connected by boundary throats. 
        pnum_dif = self._net.get_num_pores()-pnum_orig
        self._net.pore_properties['coords']=np.append(self._net.pore_properties['coords'],np.zeros((pnum_dif,3)),0)
        btlist = self._net.throat_properties['numbering'][self._net.throat_properties['type']>0]
        btnum = np.size(btlist)  
        boundary_pore_list = np.zeros((btnum,2),dtype=np.int32)
        for i in range(btnum):
            boundary_pore_list[i] = self._net.get_connected_pores(btlist[i])
        #   For each boundary pore in the pair, adopt the internal pore's coords
        for i in boundary_pore_list:
            if i[0] >= pnum_orig:
                self._net.pore_properties['coords'][i[0]] = self._net.pore_properties['coords'][i[1]]
            if i[1] >= pnum_orig:
                self._net.pore_properties['coords'][i[1]] = self._net.pore_properties['coords'][i[0]]
        #   Make lists of pores on each boundary
        face1pores = np.nonzero(self._net.pore_properties['type']==1)[0]
        face2pores = np.nonzero(self._net.pore_properties['type']==2)[0]
        face3pores = np.nonzero(self._net.pore_properties['type']==3)[0]
        face4pores = np.nonzero(self._net.pore_properties['type']==4)[0]
        face5pores = np.nonzero(self._net.pore_properties['type']==5)[0]
        face6pores = np.nonzero(self._net.pore_properties['type']==6)[0]
        #   Appropriately add or subtract a lattice constant from the appropriate
        #   dimention in the boundary pore's 'coords' value.
        for i in face1pores:
            self._net.pore_properties['coords'][i][2] += -self._Lc
        for i in face2pores:
            self._net.pore_properties['coords'][i][0] += -self._Lc
        for i in face3pores:
            self._net.pore_properties['coords'][i][1] += -self._Lc
        for i in face4pores:
            self._net.pore_properties['coords'][i][1] += self._Lc
        for i in face5pores:
            self._net.pore_properties['coords'][i][0] += self._Lc
        for i in face6pores:
            self._net.pore_properties['coords'][i][2] += self._Lc
        self._net.pore_properties['coords'][self._net.pore_properties['type']==1]
        #Update network
        self._net.update()
        
        self._logger.debug("add_boundaries: End of method")
        
    def add_opposing_boundaries(self,face,periodic=0):
        self._logger.debug("add_opposing_boundaries: Start of method")
        
        Nx = self._Nx
        Ny = self._Ny
        Nz = self._Nz
        Lx = self._Lx
        Ly = self._Ly
        Lz = self._Lz
        Lc = self._Lc
        Np = self._net.get_num_pores()
        Nt = self._net.get_num_throats()
        
        col = [-1,2,0,1,1,0,2] #Column to use in coord
        coord = [-1, Lc/2, Lc/2, Lc/2, Ly-Lc/2, Lx-Lc/2, Lz-Lc/2]
        coordperi = [-1, Lz-Lc/2, Lx-Lc/2, Ly-Lc/2, Lc/2, Lc/2, Lc/2]
        NpFace = [-1, Nx*Ny, Ny*Nz, Nx*Nz, Nx*Nz, Ny*Nz, Nx*Ny]

        #Extract pore numbers from opposing faces of the network
        tpore1 = self._net.pore_properties['numbering'][self._net.pore_properties['coords'][:,col[face]]==coord[face]]
        tpore2 = self._net.pore_properties['numbering'][self._net.pore_properties['coords'][:,col[face]]==coordperi[face]]
    
        if periodic:
            #If periodic simply link faces together
            conns = np.vstack((tpore1,tpore2)).T
            #Add elements to throat lists
            self._net.throat_properties['connections'] = np.concatenate((self._net.throat_properties['connections'],conns),axis=0)
            self._net.throat_properties['numbering'] = np.concatenate((self._net.throat_properties['numbering'],np.arange(Nt,Nt+NpFace[face],dtype=np.int32)),axis=0)
            self._net.throat_properties['type'] = np.concatenate((self._net.throat_properties['type'],np.ones(NpFace[face],dtype=np.int8)*face),axis=0)
        else:
            #If not periodic, then
            tpore1 = np.concatenate((tpore1,tpore2),axis=0)
            tpore2 = np.arange(Np,Np+2*NpFace[face])
            conns = np.vstack((tpore1,tpore2)).T
            #Add new elements to throat lists
            self._net.throat_properties['connections'] = np.concatenate((self._net.throat_properties['connections'],conns),axis=0)
            self._net.throat_properties['numbering'] = np.concatenate((self._net.throat_properties['numbering'],np.arange(Nt,Nt+2*NpFace[face],dtype=np.int32)),axis=0)
            self._net.throat_properties['type'] = np.concatenate((self._net.throat_properties['type'],np.ones(NpFace[face],dtype=np.int8)*face),axis=0)
            self._net.throat_properties['type'] = np.concatenate((self._net.throat_properties['type'],np.ones(NpFace[face],dtype=np.int8)*(7-face)),axis=0)            
            #Add new elements to pore lists
            self._net.pore_properties['numbering'] = np.concatenate((self._net.pore_properties['numbering'],np.arange(Np,Np+2*NpFace[face],dtype=np.int32)),axis=0)       
            self._net.pore_properties['type'] = np.concatenate((self._net.pore_properties['type'],np.ones(NpFace[face],dtype=np.int8)*(face)),axis=0)
            self._net.pore_properties['type'] = np.concatenate((self._net.pore_properties['type'],np.ones(NpFace[face],dtype=np.int8)*(7-face)),axis=0)
            self._net.update()
            bnum1 = self._net.pore_properties['numbering'][self._net.pore_properties['type']>(face)]
            bnum2 = self._net.pore_properties['numbering'][self._net.pore_properties['type']>(7-face)]
            pnum1 = self._net.get_neighbor_pores(bnum1,[0])
            pnum2 = self._net.get_neighbor_pores(bnum2,[0])
#            self._net.pore_properties['coords'][bnum1] = self._net.pore_properties['coords'][pnum1] - []
#            self._net.pore_properties['coords'][bnum2] = 
            
        self._logger.debug("add_opposing_boundaries: End of method")
        
      '''  
if __name__ == '__main__':
    test=Cubic(lattice_spacing=1.0,domain_size=[3,3,3],loggername='TestCubic')
    pn = test.generate()