"""
module __Cubic__: Generate simple cubic networks
==========================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import numpy as np
import scipy.stats as spst

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
        
    def generate(self,**params):
        '''
        Create Cubic network. Returns OpenPNM.Network.GenericNetwork() object.

        Parameters
        ----------

        Critical\n
        domain_size : [float,float,float]
            domain_size = [3.0,3.0,3.0] (default)\n
            Bounding cube for internal pore positions\n
        lattice_spacing : [float]
            lattice_spacing = [1.0] (default)\n
            Distance between pore centers\n
        divisions : [int,int,int]
            divisions = [3,3,3]\n
            Number of internal pores in each dimension.\n
            (Optional input. Replaces one of the above.)\n

        Optional\n
        stats_pores : dictionary
            stats_pores = {'name':'weibull_min','shape':1.5,'loc': 6e-6,'scale':2e-5} (default)\n
            Probablity distributions for random pore size assignment\n
        stats_throats : dictionary
            stats_throats = {'name':'weibull_min','shape':1.5,'loc': 6e-6,'scale':2e-5} (default)\n
            Probablity distributions for random throat size assignment\n
        btype : [logical,logical,logical]
            btype = [0,0,0] (default)\n
            Automatically create periodic throats between opposite x, y, or z faces

        Examples:
        ---------

        generate default 100x100x10 cubic network with no periodic boundaries

        >>> import OpenPNM as PNM
        >>> pn=PNM.Geometry.Cubic(domain_size=[100,100,10],lattice_spacing = 1.0)
        '''
        super(Cubic,self).generate(**params)    
        return self._net

    def _generate_setup(self,   domain_size = [],
                                divisions = [],
                                lattice_spacing = [],
                                btype = [0,0,0],
                                **params):
        r"""
        Perform applicable preliminary checks and calculations required for generation
        """
        self._logger.debug("generate_setup: Perform preliminary calculations")
        #Parse the given network size variables
        self._logger.info("Find network dimensions")
        self._btype = btype
        if domain_size and lattice_spacing and not divisions:
            self._Lc = np.float(lattice_spacing[0])
            self._Lx = np.float(domain_size[0])
            self._Ly = np.float(domain_size[1])
            self._Lz = np.float(domain_size[2])
            self._Nx = np.int(self._Lx/self._Lc)
            self._Ny = np.int(self._Ly/self._Lc)
            self._Nz = np.int(self._Lz/self._Lc)
        elif divisions and lattice_spacing and not domain_size:
            self._Lc = np.float(lattice_spacing[0])
            self._Nx = np.int(divisions[0])
            self._Ny = np.int(divisions[1])
            self._Nz = np.int(divisions[2])
            self._Lx = np.float(self._Nx*self._Lc)
            self._Ly = np.float(self._Ny*self._Lc)
            self._Lz = np.float(self._Nz*self._Lc)
        elif domain_size and divisions and not lattice_spacing:
            self._Lc = np.min(np.array(domain_size,dtype=np.float)/np.array(divisions,dtype=np.float))
            self._Nx = np.int(divisions[0])
            self._Ny = np.int(divisions[1])
            self._Nz = np.int(divisions[2])
            self._Lx = np.float(self._Nx*self._Lc)
            self._Ly = np.float(self._Ny*self._Lc)
            self._Lz = np.float(self._Nz*self._Lc)
        elif not domain_size and not divisions and not lattice_spacing:
            self._Lc = np.float(1)
            self._Lx = np.float(3)
            self._Ly = np.float(3)
            self._Lz = np.float(3)
            self._Nx = np.int(self._Lx/self._Lc)
            self._Ny = np.int(self._Ly/self._Lc)
            self._Nz = np.int(self._Lz/self._Lc)
        else:
            self._logger.error("Exactly two of domain_size, divisions and lattice_spacing must be given")
            raise Exception('Exactly two of domain_size, divisions and lattice_spacing must be given')

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
        a = np.array(np.unravel_index(ind, dims=(Nx, Ny, Nz), order='F')).T
        self._net.pore_properties['coords'] = Lc/2+Lc*np.array(np.unravel_index(ind, dims=(Nx, Ny, Nz), order='F'),dtype=np.float).T
        self._net.pore_properties['numbering'] = ind
        self._net.pore_properties['type']= np.zeros((Np,),dtype=np.int8)

        self._logger.debug("generate_pores: End of method")

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

    def _add_boundaries(self):
        r"""
        Add boundaries to network
        """
        self._logger.debug("add_boundaries: Start of method")
        #Remove all items pertaining to previously defined boundaries (if any)
        Np = self._net.get_num_pores([0])
        Nt = self._net.get_num_throats([0])
#        for item in self._net.pore_properties.keys():
#            self._net.pore_properties[item] = self._net.pore_properties[item][0:Np]
#        for item in self._net.throat_properties.keys():
#            self._net.throat_properties[item] = self._net.throat_properties[item][0:Nt]
        pnum_orig = self._net.get_num_pores([0])
        self._add_opposing_boundaries(face=2,periodic=self._btype[0]) #x faces
        self._add_opposing_boundaries(face=3,periodic=self._btype[1]) #y faces
        self._add_opposing_boundaries(face=1,periodic=self._btype[2]) #z faces

        pnum_added = self._net.get_num_pores() - pnum_orig
        self._net.pore_properties['coords'] = np.concatenate((self._net.pore_properties['coords'],np.zeros((pnum_added,3))),axis=0)
        #Add 'coords' to boundaries
        #   Generate an Nx2 array, named "boundary_pore_list" that names all
        #   pairs of pores connected by boundary throats.
        pnum_dif = self._net.get_num_pores()-pnum_orig
        btlist = self._net.throat_properties['numbering'][self._net.throat_properties['type']>0]
#        self._net.pore_properties['coords']=np.append(self._net.pore_properties['coords'],np.zeros((pnum_dif,3)),0)
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
        #Update network
        self._net.update()

        self._logger.debug("add_boundaries: End of method")

    def _add_opposing_boundaries(self,face,periodic=0):
        r"""
        Add boundaries by adding opposing faces, one pair at a time.
        """
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
        coord = np.array([-1, Lc/2, Lc/2, Lc/2, Ly-Lc/2, Lx-Lc/2, Lz-Lc/2],dtype=np.float)
        coordperi = np.array([-1, Lz-Lc/2, Lx-Lc/2, Ly-Lc/2, Lc/2, Lc/2, Lc/2],dtype=np.float)
        NpFace = np.array([-1, Nx*Ny, Ny*Nz, Nx*Nz, Nx*Nz, Ny*Nz, Nx*Ny],dtype=np.int)

        #Extract pore numbers from opposing faces of the network
        tpore1 = self._net.pore_properties['numbering'][np.abs(self._net.pore_properties['coords'][:,col[face]]-coord[face]) < np.min([Lx,Ly,Lz])/1000000]
        tpore2 = self._net.pore_properties['numbering'][np.abs(self._net.pore_properties['coords'][:,col[face]]-coordperi[face]) < np.min([Lx,Ly,Lz])/1000000]

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
#            bnum1 = self._net.pore_properties['numbering'][self._net.pore_properties['type']>(face)]
#            bnum2 = self._net.pore_properties['numbering'][self._net.pore_properties['type']>(7-face)]
#            pnum1 = self._net.get_neighbor_pores(bnum1,[0])
#            pnum2 = self._net.get_neighbor_pores(bnum2,[0])
#            self._net.pore_properties['coords'][bnum1] = self._net.pore_properties['coords'][pnum1] - []
#            self._net.pore_properties['coords'][bnum2] =

        self._logger.debug("add_opposing_boundaries: End of method")


if __name__ == '__main__':
    test=Cubic(loggername='TestCubic')
    pn = test.generate(lattice_spacing=1.0,domain_size=[3,3,3], btype=[1,1,0])