"""
module __GenericNetwork__: Abstract class to construct pore networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Network.__init__.py' file.

"""

import sys, os, collections
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM
import OpenPNM.Utilities.misc as misc
import numpy as np
import scipy as sp
import scipy.sparse as sprs
import scipy.spatial as sptl
import scipy.signal as spsg


class GenericNetwork(OpenPNM.Utilities.Tools):
    r"""
    GenericNetwork - Base class to construct pore networks

    This class contains the interface definition for the construction of networks

    Parameters
    ----------
    name : string
        Unique name for Network object
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
    loggername : string (optional)
        Define the logger name to be used on console output. Defaults to class name.

    """

    def __init__(self, name=None,**kwargs):
        r"""
        Initialize
        """
        super(GenericNetwork,self).__init__(**kwargs)
        self._logger.info("Construct Network")

        #Initialize fluid, physics, and geometry tracking lists
        self._fluids = []
        self._geometries = []
        #Initialize adjacency and incidence matrix dictionaries
        self._incidence_matrix = {}
        self._adjacency_matrix = {}
        self.reset_graphs()
        self._logger.debug("Construction of Network container")
        self.name = name
        
    def generate(self, **params):
        r"""
        Generate the network
        """
        raise NotImplementedError()
        
    #--------------------------------------------------------------------------
    '''Graph Theory and Network Query Methods'''
    #--------------------------------------------------------------------------
    def create_adjacency_matrix(self,data=None,sprsfmt='coo',dropzeros=True,sym=True):
        r"""
        Generates a weighted adjacency matrix in a desired sparse storage format

        Parameters
        ----------
        data : array_like, optional
            An array containing the throat values to enter into the matrix (In
            graph theory these are known as the 'weights').  If omitted, ones 
            are used to create a standard adjacency matrix representing 
            connectivity only.  
        sprsfmt : string, optional
            The sparse storage format to return.  Options are:
            
            * 'coo' : (default) This is the native format of OpenPNMs data
            
            * 'lil' : Enables row-wise slice of data
            
            * 'csr' : Favored by most linear algebra routines
            
        dropzeros : boolean, optional
            Remove 0 elements from the values, instead of creating 0-weighted 
            links, the default is True.
        sym : Boolean, optional
            Makes the matrix symmetric about the diagonal, the default is true.

        Returns
        -------
        Returns an adjacency matrix in the specified Scipy sparse format
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> vals = sp.rand(pn.num_throats(),) < 0.5
        >>> temp = pn.create_adjacency_matrix(data=vals,sprsfmt='csr')

        """
        self._logger.debug('create_adjacency_matrix: Start of method')
        Np   = self.num_pores()
        Nt   = self.num_throats()
        
        #Check if provided data is correct
        if data == None:
            data = sp.ones((self.num_throats(),))
        elif sp.shape(data)[0] != Nt:
            raise Exception('Received dataset of incorrect length')
            
        #Clear any zero-weighted connections
        if dropzeros:
            ind = data>0
        else:
            ind = sp.ones_like(data,dtype=bool)
            
        #Get connectivity info from network
        conn = self['throat.conns'][ind]
        row  = conn[:,0]
        col  = conn[:,1]
        data = data[ind]
        
        if sym: #Append row & col to each other, and data to itself
            row  = sp.append(row,conn[:,1])
            col  = sp.append(col,conn[:,0])
            data = sp.append(data,data)
        
        #Generate sparse adjacency matrix in 'coo' format
        temp = sprs.coo_matrix((data,(row,col)),(Np,Np))
        
        #Convert to requested format
        if sprsfmt == 'coo':
            pass #temp is already in coo format
        if sprsfmt == 'csr':
            temp = temp.tocsr()
        if sprsfmt == 'lil':
            temp = temp.tolil()
        return temp

    def create_incidence_matrix(self,data=None,sprsfmt='coo',dropzeros=True):
        r"""
        Creates an incidence matrix filled with supplied throat values

        Parameters
        ----------
        data : array_like, optional
            An array containing the throat values to enter into the matrix (In
            graph theory these are known as the 'weights').  If omitted, ones 
            are used to create a standard incidence matrix representing 
            connectivity only. 
        sprsfmt : string, optional
            The sparse storage format to return.  Options are:
            
            * 'coo' : (default) This is the native format of OpenPNMs data
            
            * 'lil' : Enables row-wise slice of data
            
            * 'csr' : Favored by most linear algebra routines
            
        dropzeros : Boolean, optional
            Remove 0 elements from values, instead of creating 0-weighted 
            links, the default is True.

        Returns
        -------
        An incidence matrix (a cousin to the adjacency matrix, useful for 
        finding throats of given a pore)

        Examples
        --------
        >>> print('nothing yet')
        nothing yet
        """
        self._logger.debug('create_incidence_matrix: Start of method')

        Nt = self.num_throats()
        Np = self.num_pores()

        #Check if provided data is valid
        if data == None:
            data = sp.ones((self.num_throats(),))
        elif sp.shape(data)[0] != Nt:
            raise Exception('Received dataset of incorrect length')

        if dropzeros:
            ind = data > 0
        else:
            ind = sp.ones_like(data, dtype=bool)

        conn = self['throat.conns'][ind]
        row  = conn[:,0]
        row = sp.append(row,conn[:,1])
        col = self.get_throat_indices('all')[ind]
        col = sp.append(col,col)
        data = sp.append(data[ind],data[ind])

        temp = sprs.coo.coo_matrix((data,(row,col)),(Np,Nt))

        #Convert to requested format
        if sprsfmt == 'coo':
            pass #temp is already in coo format
        if sprsfmt == 'csr':
            temp = temp.tocsr()
        if sprsfmt == 'lil':
            temp = temp.tolil()
        return temp
            
    def find_connected_pores(self,throats=[],flatten=False):
        r"""
        Return a list of pores connected to a list of throats

        Parameters
        ----------
        throats : array_like
            List of throats numbers
        flatten : boolean, optional
            If flatten is True (default) a 1D array of unique pore numbers
            is returned. If flatten is False each location in the the returned 
            array contains a sub-arras of neighboring pores for each input 
            throat, in the order they were sent.

        Returns
        -------
        1D array (if flatten is True) or ndarray of arrays (if flatten is False)

        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.find_connected_pores(throats=[0,1])
        array([[0, 1],
               [0, 5]])
        >>> pn.find_connected_pores(throats=[0,1],flatten=True)
        array([0, 1, 5])
        """
        Ps = self['throat.conns'][throats]
        #Ps = [sp.asarray(x) for x in Ps if x]
        if flatten:
            Ps = sp.unique(sp.hstack(Ps))
        return Ps

    def find_connecting_throat(self,P1,P2):
        r"""
        Return the throat number connecting pairs of pores

        Parameters
        ----------
        P1 , P2 : int
            The pore numbers whose throats are sought

        Returns
        -------
        Tnum : int
            Returns throat number, or empty array if pores are not connected
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.find_connecting_throat(0,1)
        array([0])
        """
        return sp.intersect1d(self.find_neighbor_throats(P1),self.find_neighbor_throats(P2))

    def find_neighbor_pores(self,pores,mode='union',flatten=True,excl_self=True):
        r"""
        Returns a list of pores neighboring the given pore(s)

        Parameters
        ----------
        pores : array_like
            ID numbers of pores whose neighbors are sought.
        flatten : boolean, optional
            If flatten is True  a 1D array of unique pore ID numbers is 
            returned. If flatten is False the returned array contains arrays
            of neighboring pores for each input pore, in the order they were 
            sent.
        excl_self : bool, optional (Default is False)
            If this is True then the input pores are not included in the 
            returned list.  This option only applies when input pores
            are in fact neighbors to each other, otherwise they are not
            part of the returned list anyway.
        mode : string, optional
            Specifies which neighbors should be returned.  The options are: 
            
            * 'union' : All neighbors of the input pores

            * 'intersection' : Only neighbors shared by all input pores 
            
            * 'not_intersection' : Only neighbors not shared by any input pores

        Returns
        -------
        neighborPs : 1D array (if flatten is True) or ndarray of ndarrays (if 
        flatten if False)

        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.find_neighbor_pores(pores=[0,2])
        array([ 1,  3,  5,  7, 25, 27])
        >>> pn.find_neighbor_pores(pores=[0,1],excl_self=True) #Find all neighbors, excluding selves
        array([ 2,  5,  6, 25, 26])
        >>> pn.find_neighbor_pores(pores=[0,2],flatten=False)
        array([array([ 1,  5, 25]), array([ 1,  3,  7, 27])], dtype=object)
        >>> pn.find_neighbor_pores(pores=[0,2],mode='intersection') #Find only common neighbors
        array([1], dtype=int64)
        >>> pn.find_neighbor_pores(pores=[0,2],mode='not_intersection') #Exclude common neighbors
        array([ 3,  5,  7, 25, 27], dtype=int64)
        >>> pn.find_neighbor_pores(pores=[0,1],mode='union') #Find all neighbors, including selves
        array([ 0,  1,  2,  5,  6, 25, 26])
        """
        pores = sp.array(pores,ndmin=1)
        try:
            neighborPs = self._adjacency_matrix['lil'].rows[[pores]]
        except:
            self._logger.info('Creating adjacency matrix, please wait')
            temp = self.create_adjacency_matrix(sprsfmt='lil')
            self._adjacency_matrix['lil'] = temp
            neighborPs = self._adjacency_matrix['lil'].rows[[pores]]
        if [sp.asarray(x) for x in neighborPs if x] == []:
            return []
        if flatten:
            #All the empty lists must be removed to maintain data type after hstack (numpy bug?)
            neighborPs = [sp.asarray(x) for x in neighborPs if x]
            neighborPs = sp.hstack(neighborPs)
            neighborPs = sp.concatenate((neighborPs,pores))
            #Remove references to input pores and duplicates
            if mode == 'not_intersection':
                neighborPs = sp.unique(sp.where(sp.bincount(neighborPs)==1)[0])
            elif mode == 'union':
                neighborPs = sp.unique(neighborPs)
            elif mode == 'intersection':
                neighborPs = sp.unique(sp.where(sp.bincount(neighborPs)>1)[0])
            if excl_self:
                neighborPs = neighborPs[~sp.in1d(neighborPs,pores)]
        else:
            for i in range(0,sp.size(pores)):
                neighborPs[i] = sp.array(neighborPs[i])
        return sp.array(neighborPs,ndmin=1)

    def find_neighbor_throats(self,pores,mode='union',flatten=True):
        r"""
        Returns a list of throats neighboring the given pore(s)

        Parameters
        ----------
        pores : array_like
            Indices of pores whose neighbors are sought
        flatten : boolean, optional
            If flatten is True (default) a 1D array of unique throat ID numbers
            is returned. If flatten is False the returned array contains arrays
            of neighboring throat ID numbers for each input pore, in the order
            they were sent.
        mode : string, optional
            Specifies which neighbors should be returned.  The options are: 
            
            * 'union' : All neighbors of the input pores

            * 'intersection' : Only neighbors shared by all input pores 
            
            * 'not_intersection' : Only neighbors not shared by any input pores

        Returns
        -------
        neighborTs : 1D array (if flatten is True) or ndarray of arrays (if
            flatten if False)

        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.find_neighbor_throats(pores=[0,1])
        array([0, 1, 2, 3, 4, 5])
        >>> pn.find_neighbor_throats(pores=[0,1],flatten=False)
        array([array([0, 1, 2]), array([0, 3, 4, 5])], dtype=object)
        """
        #Test for existance of incidence matrix
        try:
            neighborTs = self._incidence_matrix['lil'].rows[[pores]]
        except:
            self._logger.info('Creating incidence matrix, please wait')
            temp = self.create_incidence_matrix(sprsfmt='lil')
            self._incidence_matrix['lil'] = temp
            neighborTs = self._incidence_matrix['lil'].rows[[pores]]
        if [sp.asarray(x) for x in neighborTs if x] == []:
            return []
        if flatten:
            #All the empty lists must be removed to maintain data type after hstack (numpy bug?)
            neighborTs = [sp.asarray(x) for x in neighborTs if x]
            neighborTs = sp.hstack(neighborTs)
            #Remove references to input pores and duplicates
            if mode == 'not_intersection':
                neighborTs = sp.unique(sp.where(sp.bincount(neighborTs)==1)[0])
            elif mode == 'union':
                neighborTs = sp.unique(neighborTs)
            elif mode == 'intersection':
                neighborTs = sp.unique(sp.where(sp.bincount(neighborTs)>1)[0])
        else:
            for i in range(0,sp.size(pores)):
                neighborTs[i] = sp.array(neighborTs[i])
        return sp.array(neighborTs,ndmin=1)

    def num_neighbors(self,pores,flatten=False):
        r"""
        Returns an ndarray containing the number of neigbhor pores for each 
        element in pores

        Parameters
        ----------
        pores : array_like
            Pores whose neighbors are to be counted
        flatten : boolean (optional)
            If False (default) the number pore neighbors for each input are
            returned as an array.  If True the sum total number of unique 
            neighbors is counted, not including the input pores even if they 
            neighbor each other.  

        Returns
        -------
        num_neighbors : 1D array with number of neighbors in each element

        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.num_neighbors(pores=[0,1],flatten=False)
        array([3, 4], dtype=int8)
        >>> pn.num_neighbors(pores=[0,1],flatten=True)  # Sum excludes pores 0 & 1
        5
        >>> pn.num_neighbors(pores=[0,2],flatten=True)  # Sum includes pore 1, but not 0 & 2
        6
        """

        #Count number of neighbors
        if flatten:
            neighborPs = self.find_neighbor_pores(pores,flatten=True,mode='union',excl_self=True)
            num = sp.shape(neighborPs)[0]
        else:
            neighborPs = self.find_neighbor_pores(pores,flatten=False)
            num = sp.zeros(sp.shape(neighborPs),dtype=sp.int8)
            for i in range(0,sp.shape(num)[0]):
                num[i] = sp.size(neighborPs[i])
        return num
        
    def find_interface_throats(self,labels=[]):
        r'''
        Finds the throats that join two pore labels.  
        
        Parameters
        ----------
        labels : list of strings
            The labels of the two pore groups whose interface is sought
            
        Returns
        -------
        An array of throat numbers that connect the given pore groups
        
        Notes
        -----
        This method is meant to find interfaces between TWO groups, regions or 
        clusters of pores (as defined by their label).  If the input labels 
        overlap or are not adjacent, an empty array is returned. 
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.set_pore_info(label='domain1',locations=[0,1,2])
        >>> pn.set_pore_info(label='domain2',locations=[5,6,7])
        >>> pn.find_interface_throats(labels=['domain1','domain2'])
        array([1, 4, 7])
        '''
        Tind = sp.array([],ndmin=1)
        if sp.shape(labels)[0] != 2:
            self._logger.error('Exactly two labels must be given')
        else:
            P1 = self.get_pore_indices(labels=labels[0])
            P2 = self.get_pore_indices(labels=labels[1])
            #Check if labels overlap
            if sp.sum(sp.in1d(P1,P2)) > 0: 
                self._logger.error('Some labels overlap, iterface cannot be found')
            else:
                T1 = self.find_neighbor_throats(P1)
                T2 = self.find_neighbor_throats(P2)
                Tmask = sp.in1d(T1,T2)
                Tind = T1[Tmask]
        return Tind
        
    def update_network(self):
        r'''
        Regenerates the adjacency and incidence matrices
        '''
        self.create_adjacency_matrix()
        self.create_incidence_matrix()
        
    def reset_graphs(self):
        r'''
        Clears the adjacency and incidence matrices.  This is necessary after
        any manipulations functions such as trim, extend, clone, etc.
        '''
        #Re-initialize adjacency and incidence matrix dictionaries
        self._logger.debug('Resetting adjacency and incidence matrices')
        self._adjacency_matrix['coo'] = {}
        self._adjacency_matrix['csr'] = {}
        self._adjacency_matrix['lil'] = {}
        self._incidence_matrix['coo'] = {}
        self._incidence_matrix['csr'] = {}
        self._incidence_matrix['lil'] = {}
        
    #--------------------------------------------------------------------------
    '''Object Association Related Methods'''
    #--------------------------------------------------------------------------
    def regenerate_fluids(self):
        r'''
        '''
        fluids = []
        for item in self._fluids:
            fluids.append(item)
        for item in fluids:
            self._logger.info('Regenerating properties for '+item.name)
            for prop in item._prop_list:
                item.regenerate(prop_list=prop)

    def regenerate_physics(self):
        r'''
        '''
        physics = []
        for item1 in self._fluids:
            for item2 in item1._physics:
                physics.append(item2)
        for item in physics:
            self._logger.info('Regenerating pore properties for '+item.name)
            for prop in item._prop_list:
                item.regenerate(prop_list=prop)
                
    def regenerate_geometries(self):
        r'''
        '''
        #Regenerate pores first
        for item in self._geometries:
            self._logger.info('Regenerating pore properties for '+item.name)
            for prop in item._prop_list:
                if prop.split('_')[0] == 'pore':
                    item.regenerate(prop_list=prop)
        #Regenerate throats second
        for item in self._geometries:
            self._logger.info('Regenerating throat properties for '+item.name)
            for prop in item._prop_list:
                if prop.split('_')[0] == 'throat':
                    item.regenerate(prop_list=prop)
            
    def add_geometry(self,name=None,subclass='GenericGeometry',**kwargs):
        r'''
        This is a helper function for creating Geometry objects.
        
        Parameters
        ----------
        name : string, optional
            Name to apply to the newly created geometry object.  If no name is
            supplied, one will be randomly generated.
        subclass : string, optional
            The subclass to use.  If nothing is specified the GenericGeometry
            class is used.
            
        Returns
        -------
        A handle to the generated OpenPNM Geometry object.  If this handle is 
        not stored, it can be found under pn._geometries.
        
        Notes
        -----
        This method is equivalent to calling the OpenPNM module directly as:
        >>> pn = OpenPNM.Network.TestNet()
        >>> geo = OpenPNM.Geometry.GenericGeometry(network=pn,name='geom1')
        '''
        temp = OpenPNM.Geometry.__getattribute__(subclass)
        return temp(network=self,name=name,**kwargs)

    def add_fluid(self,name=None,subclass='GenericFluid',**kwargs):
        r'''
        This is a helper function for creating Fluid objects.
        
        Parameters
        ----------
        name : string, optional
            Name to apply to the newly created fluid object.  If no name is
            supplied, one will be randomly generated.
        subclass : string, optional
            The subclass to use.  If nothing is specified the GenericFluid
            class is used.
            
        Returns
        -------
        A handle to the generated OpenPNM Fluid object.  If this handle is 
        not stored, it can be found under pn._fluids.
        
        Notes
        -----
        This method is equivalent to calling the OpenPNM module directly as:
        >>> pn = OpenPNM.Network.TestNet()
        >>> fluid = OpenPNM.Fluids.GenericFluid(network=pn,name='fluid1')
        '''
        temp = OpenPNM.Fluids.__getattribute__(subclass)
        return temp(network=self,name=name,**kwargs)
        
    def add_physics(self,fluid,geometry,name=None,subclass='GenericPhysics',**kwargs):
        r'''
        '''
        temp = OpenPNM.Physics.__getattribute__(subclass)
        return temp(network=self,name=name,fluid=fluid,geometry=geometry,**kwargs)
        
    def geometries(self,name=''):
        r'''
        Retrieves Geoemtry assocaiated with the network
        
        Parameters
        ----------
        name : string, optional
            The name of the Geometry object to retrieve.  
        Returns
        -------
            If name is NOT provided, then a list of Geometry names is returned. 
            If a name IS provided, then the Geometry object of that name is 
            returned.
        '''
        if name == '':
            geoms = []
            for item in self._geometries:
                geoms.append(item.name)
        else:
            geoms = self.find_object(obj_name=name)
        return geoms
        
    def fluids(self,name=''):
        r'''
        Retrieves Fluids assocaiated with the network
        
        Parameters
        ----------
        name : string, optional
            The name of the Fluid object to retrieve.  
        Returns
        -------
            If name is NOT provided, then a list of fluid names is returned. If
            a name IS provided, then the Fluid object of that name is returned.
        '''
        if name == '':
            fluids = []
            for item in self._fluids:
                fluids.append(item.name)
        else:
            fluids = self.find_object(obj_name=name)
        return fluids
        
    #--------------------------------------------------------------------------
    '''Network Manipulation Methods'''
    #--------------------------------------------------------------------------
    def clone(self,pores,apply_label=['clone'],mode='parents'):
        r'''
        Clones the specified pores and adds them to the network
        
        Parameters
        ----------
        pores : array_like
            List of pores to clone
        apply_labels : string, or list of strings
            The labels to apply to the clones, default is 'clone'
        mode : string
            Controls the connections between parents and clones.  Options are:
            
            - 'parents': (Default) Each clone is connected only the its parent
            - 'siblings': Clones are only connected to each other in the same manner as parents were connected
            - 'isolated': No connections between parents or siblings
        '''
        if sp.shape(self.props('pore'))[0] > 1:
            self._logger.warning('Cloning an active network is messy')
        self._logger.debug(sys._getframe().f_code.co_name+': Cloning pores')
        apply_label = list(apply_label)
        #Clone pores
        Np = self.num_pores()
        Nt = self.num_throats()
        parents = sp.array(pores,ndmin=1)
        pcurrent = self.get_pore_data(prop='coords')
        pclone = pcurrent[pores,:]
        pnew = sp.concatenate((pcurrent,pclone),axis=0)
        Npnew = sp.shape(pnew)[0]
        clones = sp.arange(Np,Npnew)
        #Add clone labels to network
        for item in apply_label:
            try: self['pore.'+item]
            except: self['pore.'+item] = sp.zeros((self.num_pores(),),dtype=bool)
            try: self['throat.'+item]
            except: self['throat.'+item] = sp.zeros((self.num_throats(),),dtype=bool)
        #Add connections between parents and clones
        if mode == 'parents':
            tclone = sp.vstack((parents,clones)).T
            self.extend(pore_coords=pclone,throat_conns=tclone)
        if mode == 'siblings':
            ts = self.find_neighbor_throats(pores=pores,mode='intersection')
            tclone = self['throat.conns'][ts] + self.num_pores()
            self.extend(pore_coords=pclone,throat_conns=tclone)
        if mode == 'isolated':
            self.extend(pore_coords=pclone)
        #Apply provided labels to cloned pores
        for item in apply_label:
            self['pore.'+item][self.pores('all')>=Np] = True
            self['throat.'+item][self.throats('all')>=Nt] = True
                
        # Any existing adjacency and incidence matrices will be invalid
        self.reset_graphs()
        
    def stitch(self,heads,tails):
        r'''
        Adds throat connections between specified pores
        
        Parameters
        ----------
        pores : array_like
            An Np x 2 array contains pairs of pores that should be connected
        '''
        tc = sp.vstack((heads,tails)).T
        self.extend(throat_conns=tc)

    def extend(self,pore_coords=[],throat_conns=[]):
        r'''
        Add pores (or throats) to the network
        
        Parameters
        ----------
        pore_coords : array_like
            The coordinates of the pores to add
        throat_conns : array_like
            The throat connections to add
        
        '''
        self._logger.debug(sys._getframe().f_code.co_name+': Extending network')
        Nt = self.num_throats() + int(sp.size(throat_conns)/2)
        Np = self.num_pores() + int(sp.size(pore_coords)/3)
        #Adjust 'all' labels
        del self['pore.all'], self['throat.all']
        self['pore.all'] = sp.ones((Np,),dtype=bool)
        self['throat.all'] = sp.ones((Nt,),dtype=bool)
        #Add coords and conns
        if pore_coords != []:
            coords = sp.vstack((self['pore.coords'],pore_coords))
            self['pore.coords'] = coords
        if throat_conns != []:
            conns = sp.vstack((self['throat.conns'],throat_conns))
            self['throat.conns'] = conns
        for item in self.keys():
            if item.split('.')[1] not in ['coords','conns','all']:
                if item.split('.')[0] == 'pore':
                    N = Np
                else:
                    N = Nt
                if self[item].dtype == bool:
                    temp = self[item]
                    self[item] = sp.zeros((N,),dtype=bool)
                    self[item][temp] = True
                else:
                    temp = self[item]
                    self[item] = sp.ones((N,),dtype=float)*sp.nan
                    self[item][sp.arange(0,sp.shape(temp)[0])] = temp
        self.reset_graphs()
        
    def trim(self, pores=[], throats=[], check_health=False):
        '''
        Remove pores (or throats) from the network
        
        Parameters
        ----------
        pores (or throats) : array_like
            A boolean mask of length Np (or Nt) or a list of indices of the
            pores (or throats) to be removed.

        Notes
        -----
        The prune affects the ~selected~ pores. If you want to remove all pores
        where the diameter is less than 0.5, get a mask ie:
        pn.['pore.diameter'] < 0.5 # [True, False, ...]
        and send it over as an argument. 
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.count()
        {'pore': 125, 'throat': 300}
        >>> pn.trim(pores=[1])
        >>> pn.count()
        {'pore': 124, 'throat': 296}  # 1 pore and it's 6 throats are missing
        
        '''
        pores = np.ravel(pores)
        throats = np.ravel(throats)
        
        r'''
        TODO: This logic works but can be shorted as done in subnet
        '''
        if sp.shape(pores)[0]>0:
            Pdrop = sp.zeros((self.num_pores(),),dtype=bool)
            Pdrop[pores] = True
            Pkeep = ~Pdrop
            Tdrop = sp.zeros((self.num_throats(),),dtype=bool)
            Ts = self.find_neighbor_throats(pores)
            Tdrop[Ts] = 1
            Tkeep = ~Tdrop
        elif sp.shape(throats)[0]>0:
            Tdrop = sp.zeros((self.num_throats(),),dtype=bool)
            Tdrop[throats] = 1
            Tkeep = ~Tdrop
            Pkeep = self.get_pore_indices(labels='all')
            Pkeep = self.tomask(pores=Pkeep)
        else:
            self._logger.warning('No pores or throats recieved')
            return

        #Remap throat connections
        Pnew = sp.arange(0,sum(Pkeep),dtype=int)
        Pmap = sp.ones((self.num_pores(),),dtype=int)*-1
        Pmap[Pkeep] = Pnew
        tpore1 = self['throat.conns'][:,0]
        tpore2 = self['throat.conns'][:,1]
        Tnew1 = Pmap[tpore1[Tkeep]]
        Tnew2 = Pmap[tpore2[Tkeep]]
        
        #Adjust throat lists
        items = self.keys()
        #Write 'all' label specifically
        del self['throat.all']
        self['throat.all'] = sp.ones_like(Tnew1,dtype=bool)
        del self['pore.all']
        self['pore.all'] = sp.ones_like(Pnew,dtype=bool)
        # Write connections specifically
        del self['throat.conns']
        self['throat.conns'] = sp.vstack((Tnew1,Tnew2)).T
        # Over-write remaining data and info
        for key in items:
            if key.split('.')[1] not in ['conns','all']:
                temp = self[key]
                del self[key]
                if key.split('.')[0] == 'throat':
                    self[key] = temp[Tkeep]
                if key.split('.')[0] == 'pore':
                    self[key] = temp[Pkeep]
        
        #Reset network
        self.reset_graphs()
        
        #Check network health
        if check_health:
            self.network_health()
        
    def find_clusters(self,mask=[]):
        r'''
        Identify connected clusters of pores in the network.  
        
        Parameters
        ----------
        mask : array_like, boolean
            A list of active nodes.  This method will automatically search 
            for clusters based on site or bond connectivity depending on 
            wheather the received mask is Np or Nt long.
            
        Returns
        -------
        clusters : array_like
            An Np long list of clusters numbers
            
        '''
        if sp.shape(mask)[0] == self.num_throats():
            #Convert to boolean mask if not already
            temp = sp.zeros((self.num_throats(),),dtype=bool)
            temp[mask] = True
        elif sp.shape(mask)[0] == self.num_pores():
            conns = self.find_connected_pores(throats=self.throats())
            conns[:,0] = mask[conns[:,0]]
            conns[:,1] = mask[conns[:,1]]
            temp = sp.array(conns[:,0]*conns[:,1],dtype=bool)
        else: 
            raise Exception('Mask received was neither Nt nor Np long')
        temp = self.create_adjacency_matrix(data=temp, sprsfmt='csr', dropzeros=True)
        clusters = sprs.csgraph.connected_components(temp)[1]
        return clusters
        
    def find_duplicates(self,element='throat',mode='remove'):
        r'''
        Find duplicate pores or throats in the network
        
        Parameters
        ----------
        element : string
            Controls whethe to search for duplicate 'pore' or 'throat'
        '''
        if (self._geometries != []):
            raise Exception('Network has active Geometries, cannot proceed')
                    
        if element in ['throat','throats']:
            i = self['throat.conns'][:,0]
            j = self['throat.conns'][:,1]
            v = sp.array(self['throat.all'],dtype=int)
            Np = self.num_pores()
            temp = sprs.coo_matrix((v,(i,j)),[Np,Np])
            temp = temp.tocsr()  # Convert to CSR to combine duplicates
            temp = temp.tocoo()  # And back to COO
            if mode == 'find':
                mergedTs = sp.where(temp.data>1)
                Ps12 = sp.vstack((temp.row[mergedTs], temp.col[mergedTs])).T
                dupTs = []
                for i in range(0,sp.shape(Ps12)[0]):
                    dupTs.append(self.find_connecting_throat(Ps12[i,0],Ps12[i,1]))
                return dupTs
            elif mode == 'remove':
                self['throat.conns'] = sp.vstack((temp.row,temp.col)).T
                dict.__setitem__(self,'throat.all',sp.array(temp.data,dtype=bool))
        if element in ['pore','pores']:
            temp = misc.dist(self['pore.coords'],self['pore.coords'])
            temp = sp.triu(temp,k=1)  # Remove lower triangular of matrix
            temp = sp.where(temp==0)  # Find 0 values in distance matrix
            dupPs = temp[1]>temp[0]  # Find 0 values above diagonal
            if mode == 'find':
                return dupPs
            if mode == 'remove':
                print('not implemented yet')
                pass        
        
    def network_health(self):
        r'''
        This method check the network topological health by:
        
            (1) Checking for isolated pores
            (2) Checking for islands or isolated clusters of pores
            
        Returns
        -------
        A named tuple containing isolated pores and disconnected cluters
        '''
        #Check for individual isolated pores
        health = collections.namedtuple('network_health',['disconnected_clusters','isolated_pores'])
        health.isolated_pores = []
        health.disconnected_clusters = []
        Ps = self.num_neighbors(self.pores())
        if sp.sum(Ps==0) > 0:
            self._logger.warning(str(sp.sum(Ps==0))+' pores have no neighbors')
            health.isolated_pores = sp.where(Ps==0)[0]
        #Check for clusters of isolated pores
        Cs = self.find_clusters(self.tomask(throats=self.throats('all')))
        if sp.shape(sp.unique(Cs))[0] > 1:
            self._logger.warning('Isolated clusters exist in the network')
            for i in sp.unique(Cs):
                health.disconnected_clusters.append(sp.where(Cs==i)[0])
        return health
        
    #--------------------------------------------------------------------------
    '''Domain Geometry Methods'''
    #--------------------------------------------------------------------------
    def domain_bulk_volume(self):
        r'''
        '''
        raise NotImplementedError()

    def domain_pore_volume(self):
        r'''
        '''
        raise NotImplementedError()
        
    def domain_length(self,face_1,face_2):
        r'''
        Calculate the distance between two faces of the network
        
        Parameters
        ----------
        face_1 and face_2 : array_like
            Lists of pores belonging to opposite faces of the network
            
        Returns
        -------
        The length of the domain in the specified direction
        
        Notes
        -----
        - Does not yet check if input faces are perpendicular to each other
        '''
        #Ensure given points are coplanar before proceeding
        if misc.iscoplanar(self['pore.coords'][face_1]) and misc.iscoplanar(self['pore.coords'][face_2]):
            #Find distance between given faces
            x = self['pore.coords'][face_1]
            y = self['pore.coords'][face_2]
            Ds = misc.dist(x,y)
            L = sp.median(sp.amin(Ds,axis=0))
        else:
            self._logger.warning('The supplied pores are not coplanar. Length will be approximate.')
            f1 = self['pore.coords'][face_1]
            f2 = self['pore.coords'][face_2]
            distavg = [0,0,0]
            distavg[0] = sp.absolute(sp.average(f1[:,0]) - sp.average(f2[:,0]))
            distavg[1] = sp.absolute(sp.average(f1[:,1]) - sp.average(f2[:,1]))
            distavg[2] = sp.absolute(sp.average(f1[:,2]) - sp.average(f2[:,2]))
            L = max(distavg)
        return L
        
    def domain_area(self,face):
        r'''
        Calculate the area of a given network face
        
        Parameters
        ----------
        face : array_like
            List of pores of pore defining the face of interest
            
        Returns
        -------
        The area of the specified face
        '''
        if misc.iscoplanar(self['pore.coords'][face]):
            #Find area of inlet face
            x = self['pore.coords'][face]
            y = x
            As = misc.dist(x,y)
            temp = sp.amax(As,axis=0)
            h = sp.amax(temp)
            corner1 = sp.where(temp==h)[0][0]
            p = spsg.argrelextrema(As[corner1,:],sp.greater)[0]
            a = sp.amin(As[corner1,p])
            o = h*sp.sin(sp.arccos((a/h)))
            A = o*a
        else:
            self._logger.warning('The supplied pores are not coplanar. Area will be approximate')
            x = self['pore.coords'][face]
            y = x
            As = misc.dist(x,y)
            temp = sp.amax(As,axis=0)
            h = sp.amax(temp)
            corner1 = sp.where(temp==h)[0][0]
            p = spsg.argrelextrema(As[corner1,:],sp.greater)[0]
            a = sp.amin(As[corner1,p])
            o = h*sp.sin(sp.arccos((a/h)))
            A = o*a
        return A        
        
    #--------------------------------------------------------------------------
    '''Miscillaneous Methods'''
    #--------------------------------------------------------------------------
    def subset(self,pores,name=None,incl_labels=True,incl_props=True):
        r'''
        Create a new sub-network from a list of pores.
        
        Parameters
        ----------
        pores : array_like
            A list of pores from which to create the new network
        incl_labels : bool, default is True
            Specifies whether to keep labels from main network
        incl_props : bool, default is True
            Specifies whather to keep properties from main network
        name : string, optional
            The name to apply to the new network object
            
        Returns
        -------
        newpnm : OpenPNM Object
            Returns a new network object
            
        Notes
        -----
        This method adds a pore and throat label to the master network (pn) 
        indicating which pores and throats were part of the sub-network (sn).  
        This means that the results of the sub-network can be applied to the 
        correct pores and throats in the main network by calling 
        pn.pores(sn.name), and pn.throats(sn.name).
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.count()
        {'pore': 125, 'throat': 300}
        >>> Ps = pn.pores(['top','bottom','front'])
        >>> psn = pn.subset(pores=Ps)
        >>> sn.count()
        {'pore': 65, 'throat': 112}
        >>> sn.labels()[0:3]  # It automatically transfers labels and props
        ['pore.all', 'pore.back', 'pore.bottom']
        >>>  sn = pn.subset(pores=Ps,incl_labels=False)
        >>> sn.labels()
        ['pore.all', 'throat.all']
        '''
        newpnm = OpenPNM.Network.GenericNetwork(name=name)
        pores = sp.array(pores,ndmin=1)
        throats = self.find_neighbor_throats(pores=pores,mode='intersection',flatten=True)
        
        #Remap throats on to new pore numbering
        Pmap = sp.zeros_like(self.pores())*sp.nan
        Pmap[pores] = sp.arange(0,sp.shape(pores)[0])
        tpore1 = sp.array(Pmap[self['throat.conns'][throats,0]],dtype=int)
        tpore2 = sp.array(Pmap[self['throat.conns'][throats,1]],dtype=int)
        newpnm['throat.conns'] = sp.vstack((tpore1,tpore2)).T
        newpnm['pore.coords'] = self['pore.coords'][pores]
        
        #Now scan through labels and props, and keep if needed
        newpnm['pore.all'] = self['pore.all'][pores]
        newpnm['throat.all'] = self['throat.all'][throats]
        if incl_labels == True:
            labels = self.labels()
            labels.remove('pore.all')
            labels.remove('throat.all')
            for item in labels:
                if item.split('.')[0] == 'pore':
                    newpnm[item] = self[item][pores]
                if item.split('.')[0] == 'throat':
                    newpnm[item] = self[item][throats]
        if incl_props == True:
            props = self.props()
            props.remove('throat.conns')
            props.remove('pore.coords')
            for item in props:
                if item.split('.')[0] == 'pore':
                    newpnm[item] = self[item][pores]
                if item.split('.')[0] == 'throat':
                    newpnm[item] = self[item][throats]
        
        #Append pore and throat mapping to main network as attributes and data
        self['pore.'+newpnm.name] = sp.zeros_like(self['pore.all'],dtype=bool)
        self['pore.'+newpnm.name][pores] = True
        self['throat.'+newpnm.name] = sp.zeros_like(self['throat.all'],dtype=bool)
        self['throat.'+newpnm.name][throats] = True
        
        import types
        newpnm.subset_fluid = types.MethodType(subset_fluid, newpnm)     
        
        return newpnm
        
#------------------------------------------------------------------------------
'''Additional Functions'''
#------------------------------------------------------------------------------
# These functions are not automatically attached to the network, but can be
# using object.method_name = types.MethodType(method_name, object)
  
def subset_fluid(self,fluid):
    r'''
    This method is appended to subset networks. It takes a fluid from the main
    network, and converts it to the size and shape of the sub-network.
    
    Parameters
    ----------
    fluid : OpenPNM Fluid Object
        A fluid object that is associated with the main network from which the
        subnetwork was extracted.
        
    Returns
    -------
    newfluid : OpenPNM Fluid Object
        A fluid object with the same shape as the sub-network.  It contains all
        the data of the main fluid, but not the property calculation methods.
    '''
    newfluid = OpenPNM.Fluids.GenericFluid(network=self)
    for item in fluid.props(mode='vectors'):
        if item.split('.')[0] == 'pore':
            newfluid[item] = fluid[item][fluid._net['pore.'+self.name]]
        if item.split('.')[0] == 'throat':
            newfluid[item] = fluid[item][fluid._net['throat.'+self.name]]
    for item in fluid.props(mode='scalars'):
        if item.split('.')[0] == 'pore':
            newfluid[item] = fluid[item]
        if item.split('.')[0] == 'throat':
            newfluid[item] = fluid[item]
    return newfluid

if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)
    
    
    
    
    
    
    
    
    
    
