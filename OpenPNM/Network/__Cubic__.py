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
import sys
from OpenPNM.Network import GenericNetwork

class Cubic(GenericNetwork):
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
    >>> pn = OpenPNM.Network.Cubic(name='test_cubic').generate(lattice_spacing=[1],divisions=[5,5,5])

    """

    def __init__(self, **kwargs):

        super(Cubic,self).__init__(**kwargs)
        self._logger.debug(self.__class__.__name__+": Execute constructor")

    def generate(self,**params):
        '''
        Invokes the creation of a Cubic network

        Parameters
        ----------
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

        '''
        self._logger.info(sys._getframe().f_code.co_name+": Start of network topology generation")
        self._generate_setup(**params)
        self._generate_pores()
        self._generate_throats()
        #self.add_boundaries()
        self._add_labels()
        self._logger.debug(sys._getframe().f_code.co_name+": Network generation complete")
        return self

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
        self._logger.info(sys._getframe().f_code.co_name+": Creating specified number of pores")
        Nx = self._Nx
        Ny = self._Ny
        Nz = self._Nz
        Lc = self._Lc
        Np = Nx*Ny*Nz
        ind = np.arange(0,Np)
        pore_coords = Lc/2+Lc*np.array(np.unravel_index(ind, dims=(Nx, Ny, Nz), order='F'),dtype=np.float).T
        self.set_pore_data(prop='coords',data=pore_coords)
        self.set_pore_data(prop='numbering',data=ind)
        self.set_pore_info(prop='all',locations=np.ones_like(ind))
        self._logger.debug(sys._getframe().f_code.co_name+": End of pore creation")

    def _generate_throats(self):
        r"""
        Generate the throats (connections, numbering and types)
        """
        self._logger.info(sys._getframe().f_code.co_name+": Creating specified number of throats")

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
        self.set_throat_data(prop='connections',data=connections)      
        self.set_throat_data(prop='numbering',data=np.arange(0,np.shape(tpore1)[0]))
        self.set_throat_info(prop='all',locations=np.ones_like(np.arange(0,np.shape(tpore1)[0])))        
        self._logger.debug(sys._getframe().f_code.co_name+": End of throat creation")
        
    def _add_labels(self):
        self._logger.info(sys._getframe().f_code.co_name+": Applying labels")
        coords = self.get_pore_data(prop='coords')
        self.set_pore_info(prop='front',locations=coords[:,0]<=self._Lc)
        self.set_pore_info(prop='left',locations=coords[:,1]<=self._Lc)
        self.set_pore_info(prop='bottom',locations=coords[:,2]<=self._Lc)
        self.set_pore_info(prop='back',locations=coords[:,0]>=(self._Lc*(self._Nx-1)))
        self.set_pore_info(prop='right',locations=coords[:,1]>=(self._Lc*(self._Ny-1)))
        self.set_pore_info(prop='top',locations=coords[:,2]>=(self._Lc*(self._Nz-1)))
        self.set_pore_info(prop='internal',locations=self.get_pore_indices(),is_indices=True)
        #Add throat labels based IF both throat's neighbors have label in common
        for item in ['top','bottom','left','right','front','back']:
            ps = self.get_pore_indices(item)
            ts = self.get_neighbor_throats(ps)
            ps = self.get_connected_pores(ts)
            ps0 = self.get_pore_info(prop=item)[ps[:,0]]
            ps1 = self.get_pore_info(prop=item)[ps[:,1]]
            ts = ts[ps1*ps0]
            self.set_throat_info(prop=item,locations=ts,is_indices=True)
        self.set_throat_info(prop='internal',locations=self.get_throat_indices(),is_indices=True)
        self._logger.debug(sys._getframe().f_code.co_name+": End")

    def _generate_boundaries(self,net,**params):

        self._logger.info("generate_boundaries: Define edge points of the pore network and stitch them on")
        self._generate_setup(**params)
        Lc = self._Lc
        #Takes the original dimensions of the first network.
        params['divisions'] = [self._Nx, self._Ny, self._Nz]
        temp_divisions = sp.copy(params['divisions']).tolist()
        pore_coords = self._net.get_pore_data(prop='coords')
        for i in range(0,6):
            params['divisions'] = temp_divisions
            temp_divisions = sp.copy(params['divisions']).tolist()
            params['divisions'][i%3] = 1
            edge_net = super(Cubic, self).generate(**params) # Generate networks based on the size of the edges. These are computed in divisions.
            displacement = [0,0,0]
            if i < 3:
                displacement[i%3] = pore_coords[:,i%3].max() + 0.5*Lc # This network is properly displaced and then translated.
            else:
                displacement[i%3] = -1*(pore_coords[:,i%3].min() + 0.5*Lc)
            self.translate_coordinates(edge_net,displacement)
            self.stitch_network(net,edge_net,edge = i+1,stitch_nets = 0) # We stitch the networks to generate the new throats and append all properties.

        self._stitch_throats(net,**params)
        self._net.set_pore_data(prop='type',data=self._add_boundary_pore_type(net))
        self._net.set_throat_data(prop='type',data=self._add_boundary_throat_type(net))
        self._logger.debug("generate_boundaries: End of method")
        return net

    def stitch_network(self,net1,net2,edge = 0, stitch_nets = 1, stitch_side = [],**params):
        r"""
        Stitch two networks together

        Parameters
        ----------
        net1 : OpenPNM Network Object
            The network that is stiched to

        net2 : OpenPNM Network Object
            The network that is stitched

        edge : default value of 0. This changed if you are adding boundaries to an existing network.

        stitch_nets : default value of 1. This assumes you are stitching 2 separately generated networks together. If the value
            is changed to a value of 0, you don't append the throats of the stitched network because you are adding boundary pores.

        stitch_side : string (optional)
            When given, this flag moves net2 into position automatically

        """
        net1_pore_coords = net1.get_pore_data(prop='coords')
        net2_pore_coords = net2.get_pore_data(prop='coords')
        net1_pore_numbering = net1.get_pore_data(prop='numbering')
        net2_pore_numbering = net2.get_pore_data(prop='numbering')
        net1_pore_volume = net1.get_pore_data(prop='volume')
        net2_pore_volume = net2.get_pore_data(prop='volume')
        net1_pore_seed = net1.get_pore_data(prop='seed')
        net2_pore_seed = net2.get_pore_data(prop='seed')
        net1_pore_diameter = net1.get_pore_data(prop='diameter')
        net2_pore_diameter = net2.get_pore_data(prop='diameter')
        net1_pore_type = net1.get_pore_data(prop='type')
        net2_pore_type = net2.get_pore_data(prop='type')
        net1_throat_numbering = net1.get_throat_data(prop='numbering')
        net2_throat_numbering = net2.get_throat_data(prop='numbering')
        net1_throat_seed = net1.get_throat_data(prop='seed')
        net2_throat_seed = net2.get_throat_data(prop='seed')
        net1_throat_diameter = net1.get_throat_data(prop='diameter')
        net2_throat_diameter = net2.get_throat_data(prop='diameter')
        net1_throat_volume = net1.get_throat_data(prop='volume')
        net2_throat_volume = net2.get_throat_data(prop='volume')
        net1_throat_length = net1.get_throat_data(prop='length')
        net2_throat_length = net2.get_throat_data(prop='length')
        net1_throat_type = net1.get_throat_data(prop='type')
        net2_throat_type = net2.get_throat_data(prop='type')
        
        if stitch_nets: # This is the side where we stitch the new network onto.
            #if ~stitch_size:
                #stitch in place
            if stitch_side == 'left':
                def_translate = net1_pore_coords[:,0].max()+net1_pore_coords[:,0].min()
                OpenPNM.Geometry.GenericGeometry.translate_coordinates(net2,displacement=[def_translate,0,0])
            elif stitch_side == 'right':
                def_translate = net1_pore_coords[:,0].max()+net1_pore_coords[:,0].min()
                OpenPNM.Geometry.GenericGeometry.translate_coordinates(net2,displacement=[-1*def_translate,0,0])
            elif stitch_side == 'back':
                def_translate = net1_pore_coords[:,1].max()+net1_pore_coords[:,1].min()
                OpenPNM.Geometry.GenericGeometry.translate_coordinates(net2,displacement=[0,def_translate,0])
            elif stitch_side == 'front':
                def_translate = net1_pore_coords[:,1].max()+net1_pore_coords[:,1].min()
                OpenPNM.Geometry.GenericGeometry.translate_coordinates(net2,displacement=[0,-1*def_translate,0])
            elif stitch_side == 'top':
                def_translate = net1_pore_coords[:,2].max()+net1_pore_coords[:,2].min()
                OpenPNM.Geometry.GenericGeometry.translate_coordinates(net2,displacement=[0,0,def_translate])
            elif stitch_side == 'bottom':
                def_translate = net1_pore_coords[:,2].max()+net1_pore_coords[:,2].min()
                OpenPNM.Geometry.GenericGeometry.translate_coordinates(net2,displacement=[0,0,-1*def_translate])

        #Concatenate all conserved properties.
        net2_pore_numbering   = len(net1_pore_numbering) + net2_pore_numbering
        net1_pore_numbering   = sp.concatenate((net1_pore_numbering,net2_pore_numbering),axis=0)
        net1_pore_volume      = sp.concatenate((net1_pore_volume,net2_pore_volume),axis = 0)
        net1_pore_seed        = sp.concatenate((net1_pore_seed,net2_pore_seed),axis = 0)
        net1_pore_diameter    = sp.concatenate((net1_pore_diameter,net2_pore_diameter),axis = 0)
        net1_pore_coords      = sp.concatenate((net1_pore_coords,net2_pore_coords),axis = 0)
        net1_pore_type        = sp.concatenate((net1_pore_type,net2_pore_type),axis = 0)
        net2_throat_numbering = len(net1_throat_numbering) + net2_throat_numbering
        net1_throat_numbering = sp.concatenate((net1_throat_numbering,net2_throat_numbering),axis=0)
        net1_throat_seed      = sp.concatenate((net1_throat_seed,net2_throat_seed),axis=0)
        net1_throat_diameter  = sp.concatenate((net1_throat_diameter,net2_throat_diameter),axis=0)
        net1_throat_volume    = sp.concatenate((net1_throat_volume,net2_throat_volume),axis=0)
        net1_throat_length    = sp.concatenate((net1_throat_length,net2_throat_length),axis=0)

        if stitch_nets:
            self._stitch_throats(net1,**params)
            net1_pore_type = self._add_boundary_pore_type(net1)
            net1_throat_type = self._add_boundary_throat_type(net1)

        net1.set_pore_data(prop='coords',data=net1_pore_coords)
        net2.set_pore_data(prop='coords',data=net2_pore_coords)
        net1.set_pore_data(prop='numbering',data=net1_pore_numbering)
        net2.set_pore_data(prop='numbering',data=net2_pore_numbering)
        net1.set_pore_data(prop='volume',data=net1_pore_volume)
        net2.set_pore_data(prop='volume',data=net2_pore_volume)
        net1.set_pore_data(prop='seed',data=net1_pore_seed)
        net2.set_pore_data(prop='seed',data=net2_pore_seed)
        net1.set_pore_data(prop='diameter',data=net1_pore_diameter)
        net2.set_pore_data(prop='diameter',data=net2_pore_diameter)
        net1.set_pore_data(prop='type',data=net1_pore_type)
        net2.set_pore_data(prop='type',data=net2_pore_type)
        net1.set_throat_data(prop='numbering',data=net1_throat_numbering)
        net2.set_throat_data(prop='numbering',data=net2_throat_numbering)
        net1.set_throat_data(prop='seed',data=net1_throat_seed)
        net2.set_throat_data(prop='seed',data=net2_throat_seed)
        net1.set_throat_data(prop='diameter',data=net1_throat_diameter)
        net2.set_throat_data(prop='diameter',data=net2_throat_diameter)
        net1.set_throat_data(prop='volume',data=net1_throat_volume)
        net2.set_throat_data(prop='volume',data=net2_throat_volume)
        net1.set_throat_data(prop='length',data=net1_throat_length)
        net2.set_throat_data(prop='length',data=net2_throat_length)
        net1.set_throat_data(prop='type',data=net1_throat_type)
        net2.set_throat_data(prop='type',data=net2_throat_type)

    def _stitch_throats(self,net,
                     stats_throats = {'name' : 'weibull_min',
                                     'shape' : 1.5,
                                       'loc' : 6e-6,
                                     'scale' : 2e-5},
                          **params):
        r"""
        Stitch two networks together OR adds the boundary throats to an existing network

        Parameters
        ----------
        net : OpenPNM Network Object
            The network that is stiched, whos throats are being added.

        """
        pts = net.get_pore_data(prop='coords')
        tri = sptl.Delaunay(pts)
        I = []
        J = []
        V = []

        dist_comb = list(itr.combinations(range(4),2)) # [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]. This compares each colony with its neighbor with no repeats.

        for i in range(len(tri.simplices)):
            for j in range(len(dist_comb)):
                point_1 = tri.simplices[i,dist_comb[j][0]]
                point_2 = tri.simplices[i,dist_comb[j][1]]
                coords_1 = tri.points[point_1]
                coords_2 = tri.points[point_2]
                dist =  np.sqrt((coords_2[0] - coords_1[0]) ** 2 + (coords_2[1] - coords_1[1]) ** 2 + (coords_2[2] - coords_1[2]) ** 2)
                V.append(dist)
                I.append(point_1)
                J.append(point_2)

                    # Takes all the IJV values and puts it into a sparse matrix.
        spar_connections = sp.sparse.coo_matrix((np.array(V),(np.array(I),np.array(J))))
        ind = np.where(spar_connections.data < min(spar_connections.data) + 0.001)
        prelim_connections = np.array((spar_connections.row[ind],spar_connections.col[ind]))
        connections = np.zeros((len(prelim_connections[0]),2),np.int)

        for i in range(len(prelim_connections[0])):
            connections[i,0] = prelim_connections[0][i]
            connections[i,1] = prelim_connections[1][i]

        b = np.ascontiguousarray(connections).view(np.dtype((np.void, connections.dtype.itemsize * connections.shape[1])))
        _, idx = np.unique(b, return_index=True) # TAKEN FROM STACK OVERFLOW FOR THE TIME BEING.
        connections = connections[idx]

        net.set_throat_data(prop='connections',data=connections)        
        net.set_throat_data(prop='numbering',data=np.arange(0,len(connections[:,0])))        
        net.set_throat_data(prop='type',data=np.zeros(len(connections[:,0]),dtype=np.int8))        

        start = len(net.get_throat_data(prop='diameter'))
        end = len(net.get_throat_data(prop='connections'))
        old_seeds = net.get_throat_data(prop='seed').copy()
        new_seeds = sp.amin(net.get_pore_data(prop='seed')[net.get_throat_data(prop='connections')[start:end]],1)
        net.set_throat_data(prop='seed',data = np.concatenate((old_seeds,new_seeds)))

        prob_fn = getattr(spst,stats_throats['name'])
        P = prob_fn(stats_throats['shape'],loc=stats_throats['loc'],scale=stats_throats['scale'])
        net.set_throat_data(prop='diameter',data= P.ppf(net.get_throat_data(prop='seed')))
        net.set_throat_data(prop='length',data = sp.zeros_like(net.get_throat_data(prop='type')))
        C1 = net.get_pore_data(prop='coords')[net.get_throat_data(prop='connections')[:,0]]
        C2 = net.get_pore_data(prop='coords')[net.get_throat_data(prop='connections')[:,1]]
        E = sp.sqrt(sp.sum((C1-C2)**2,axis=1))  #Euclidean distance between pores
        D1 = net.get_pore_data(prop='diameter')[net.get_throat_data(prop='connections')[:,0]]
        D2 = net.get_pore_data(prop='diameter')[net.get_throat_data(prop='connections')[:,1]]
        net.set_throat_data(prop='length',data=E - (D1 + D2)/2)
        net.set_throat_data(prop='volume',data=net.get_throat_data(prop='length')*net.get_throat_data(prop='diameter')**2)

    def _add_boundary_throat_type(self,net):
        throat_type = np.zeros(len(net.get_throat_data(prop='type')))

        for i in range(0,len(throat_type)):
            temp1 = net.get_pore_data(prop='type')[net.get_throat_data(prop='connections')[i,0]]
            temp2 = net.get_pore_data(prop='type')[net.get_throat_data(prop='connections')[i,1]]
            if max(temp1,temp2) > 0:
                throat_type[i] = max(temp1,temp2)
        return throat_type

    def _add_boundary_pore_type(self,net):
        pore_type = np.zeros(len(net.get_pore_data(prop='type')))
        for i in range(3):
            bound_1 = net.get_pore_data(prop='coords')[:,i].min()
            bound_2 = net.get_pore_data(prop='coords')[:,i].max()
            bound_ind_1 = np.where(net.get_pore_data(prop='coords')[:,i] == bound_1)
            bound_ind_2 = np.where(net.get_pore_data(prop='coords')[:,i] == bound_2)
            pore_type[bound_ind_1] = i+1
            pore_type[bound_ind_2] = 6-i
        return pore_type

    def _add_boundaries(self):
        r"""
        Add boundaries to network
        """
#        self._logger.debug("add_boundaries: Start of method")
        #Remove all items pertaining to previously defined boundaries (if any)
#        for item in self._net.pore_data.keys():
#            self._net.pore_data[item] = self._net.pore_data[item][0:Np]
#        for item in self._net.throat_data.keys():
#            self._net.throat_data[item] = self._net.throat_data[item][0:Nt]
        pnum_orig = self._net.get_num_pores([0])
        self._add_opposing_boundaries(face=2,periodic=self._btype[0]) #x faces
        self._add_opposing_boundaries(face=3,periodic=self._btype[1]) #y faces
        self._add_opposing_boundaries(face=1,periodic=self._btype[2]) #z faces

        pnum_added = self._net.get_num_pores() - pnum_orig
        self._net.set_pore_data(prop='coords',data=np.concatenate((self._net.get_pore_data(prop='coords'),np.zeros((pnum_added,3))),axis=0))
        #Add 'coords' to boundaries
        #   Generate an Nx2 array, named "boundary_pore_list" that names all
        #   pairs of pores connected by boundary throats.
        pnum_dif = self._net.get_num_pores()-pnum_orig
        btlist = self._net.get_throat_data(prop='numbering')[self._net.get_throat_data(prop='type')>0]
#        self._net.pore_data['coords']=np.append(self._net.pore_data['coords'],np.zeros((pnum_dif,3)),0)
        btnum = np.size(btlist)
        boundary_pore_list = np.zeros((btnum,2),dtype=np.int32)
        for i in range(btnum):
            boundary_pore_list[i] = self._net.get_connected_pores(btlist[i])
        #   For each boundary pore in the pair, adopt the internal pore's coords
        for i in boundary_pore_list:
            if i[0] >= pnum_orig:
                pore_coords = self._net.get_pore_data(prop='coords')
                pore_coords[i[0]] = pore_coords[i[1]]
                self._net.set_pore_data(prop='coords',data=pore_coords)
            if i[1] >= pnum_orig:
                pore_coords = self._net.get_pore_data(prop='coords')
                pore_coords[i[1]] = pore_coords[i[0]]
                self._net.set_pore_data(prop='coords',data=pore_coords)
        #   Make lists of pores on each boundary
        face1pores = np.nonzero(self._net.get_pore_data(prop='type')==1)[0]
        face2pores = np.nonzero(self._net.get_pore_data(prop='type')==2)[0]
        face3pores = np.nonzero(self._net.get_pore_data(prop='type')==3)[0]
        face4pores = np.nonzero(self._net.get_pore_data(prop='type')==4)[0]
        face5pores = np.nonzero(self._net.get_pore_data(prop='type')==5)[0]
        face6pores = np.nonzero(self._net.get_pore_data(prop='type')==6)[0]
        #   Appropriately add or subtract a lattice constant from the appropriate
        #   dimention in the boundary pore's 'coords' value.
        for i in face1pores:
            self._net.pore_data['coords'][i][2] += -self._Lc
        for i in face2pores:
            self._net.pore_data['coords'][i][0] += -self._Lc
        for i in face3pores:
            self._net.pore_data['coords'][i][1] += -self._Lc
        for i in face4pores:
            self._net.pore_data['coords'][i][1] += self._Lc
        for i in face5pores:
            self._net.pore_data['coords'][i][0] += self._Lc
        for i in face6pores:
            self._net.pore_data['coords'][i][2] += self._Lc
        #Update network
        self._net.update()

#        self._logger.debug("add_boundaries: End of method")

    def _add_opposing_boundaries(self,face,periodic=0):
        r"""
        Add boundaries by adding opposing faces, one pair at a time.
        """
#        self._logger.debug("add_opposing_boundaries: Start of method")

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
        tpore1 = self._net.pore_data['numbering'][np.abs(self._net.pore_data['coords'][:,col[face]]-coord[face]) < np.min([Lx,Ly,Lz])/1000000]
        tpore2 = self._net.pore_data['numbering'][np.abs(self._net.pore_data['coords'][:,col[face]]-coordperi[face]) < np.min([Lx,Ly,Lz])/1000000]

        if periodic:
            #If periodic simply link faces together
            conns = np.vstack((tpore1,tpore2)).T
            #Add elements to throat lists
            self._net.throat_data['connections'] = np.concatenate((self._net.throat_data['connections'],conns),axis=0)
            self._net.throat_data['numbering'] = np.concatenate((self._net.throat_data['numbering'],np.arange(Nt,Nt+NpFace[face],dtype=np.int32)),axis=0)
            self._net.throat_data['type'] = np.concatenate((self._net.throat_data['type'],np.ones(NpFace[face],dtype=np.int8)*face),axis=0)
        else:
            #If not periodic, then
            tpore1 = np.concatenate((tpore1,tpore2),axis=0)
            tpore2 = np.arange(Np,Np+2*NpFace[face])
            conns = np.vstack((tpore1,tpore2)).T
            #Add new elements to throat lists
            self._net.throat_data['connections'] = np.concatenate((self._net.throat_data['connections'],conns),axis=0)
            self._net.throat_data['numbering'] = np.concatenate((self._net.throat_data['numbering'],np.arange(Nt,Nt+2*NpFace[face],dtype=np.int32)),axis=0)
            self._net.throat_data['type'] = np.concatenate((self._net.throat_data['type'],np.ones(NpFace[face],dtype=np.int8)*face),axis=0)
            self._net.throat_data['type'] = np.concatenate((self._net.throat_data['type'],np.ones(NpFace[face],dtype=np.int8)*(7-face)),axis=0)
            #Add new elements to pore lists
            self._net.pore_data['numbering'] = np.concatenate((self._net.pore_data['numbering'],np.arange(Np,Np+2*NpFace[face],dtype=np.int32)),axis=0)
            self._net.pore_data['type'] = np.concatenate((self._net.pore_data['type'],np.ones(NpFace[face],dtype=np.int8)*(face)),axis=0)
            self._net.pore_data['type'] = np.concatenate((self._net.pore_data['type'],np.ones(NpFace[face],dtype=np.int8)*(7-face)),axis=0)
            self._net.update()


if __name__ == '__main__':
    pn = OpenPNM.Network.Cubic(name='cubic_1',loglevel=10).generate(lattice_spacing=[1.0],domain_size=[3,3,3])
    print(pn.name)
