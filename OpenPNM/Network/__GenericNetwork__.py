"""
module __GenericNetwork__: Abstract class to construct pore networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Network.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM
import numpy as np
import scipy as sp
import scipy.sparse as sprs

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
        self.name = name
        #Initialize fluid, physics, and geometry tracking lists
        self._fluids = {}
        self._geometries = {}
        #Initialize adjacency and incidence matrix dictionaries
        self._incidence_matrix = {}
        self._adjacency_matrix = {}
        self.reset_graphs()
        self._logger.debug("Construction of Network container")

        
    def generate(self, **params):
        r"""
        Generate the network
        """
        self._logger.error('This method is not implemented')
 
    #--------------------------------------------------------------------------
    '''pore_data and throat_data interpolation methods'''
    #--------------------------------------------------------------------------
    def interpolate_pore_data(self,Tvals=None):
        r"""
        Determines a pore property as the average of it's neighboring throats

        Parameters
        ----------
        Tvals : array_like
            The array of throat information to be interpolated
        
        Returns
        -------
        Pvals : array_like
            An array of size Np contain interpolated throat data
            
        See Also
        --------
        interpolate_throat_data

        Notes
        -----
        This uses an unweighted average, without attempting to account for distances or sizes of pores and throats.
        
        Examples
        --------
        

        """
        if sp.size(Tvals)==1:
            Pvals = Tvals
        elif sp.size(Tvals) != self.num_throats():
            raise Exception('The list of throat information received was the wrong length')
        else:
            Pvals = sp.zeros((self.num_pores()))
            #Only interpolate conditions for internal pores, type=0
            Pnums = self.get_pore_indices()
            nTs = self.find_neighbor_throats(Pnums,flatten=False)
            for i in sp.r_[0:sp.shape(nTs)[0]]:
                Pvals[i] = sp.mean(Tvals[nTs[i]])
        return Pvals

    def interpolate_throat_data(self,Pvals=None):
        r"""
        Determines a throat property as the average of it's neighboring pores

        Parameters
        ----------
        Pvals : array_like
            The array of the pore condition to be interpolated
        
        Returns
        -------
        Tvals : array_like
            An array of size Nt contain interpolated pore data
            
        See Also
        --------
        interpolate_throat_data

        Notes
        -----
        This uses an unweighted average, without attempting to account for distances or sizes of pores and throats.

        Examples
        --------

        """
        if sp.size(Pvals)==1:
            Tvals = Pvals
        elif sp.size(Pvals) != self.num_pores():
            raise Exception('The list of pore information received was the wrong length')
        else:
            Tvals = sp.zeros((self.num_throats()))
            #Interpolate values for all throats, including those leading to boundary pores
            Tnums = sp.r_[0:self.num_throats()]
            nPs = self.find_connected_pores(Tnums,flatten=False)
            for i in sp.r_[0:sp.shape(nPs)[0]]:
                Tvals[i] = sp.mean(Pvals[nPs[i]])
        return Tvals

    def amalgamate_pore_data(self,fluids='all'):
        r"""
        Returns a dictionary containing ALL pore data from all fluids, physics and geometry objects
        """
        self._pore_data_amalgamate = {}
        if type(fluids)!= sp.ndarray and fluids=='all':
            fluids = self._fluids
        elif type(fluids)!= sp.ndarray: 
            fluids = sp.array(fluids,ndmin=1)
        #Add fluid data
        for item in fluids:
            if type(item)==sp.str_: item =  self.find_object_by_name(item)
            for key in item._pore_data.keys():
                if sp.amax(item._pore_data[key]) < sp.inf:
                    dict_name = item.name+'_pore_'+key
                    self._pore_data_amalgamate.update({dict_name : item._pore_data[key]})
            for key in item._pore_info.keys():
                if sp.amax(item._pore_info[key]) < sp.inf:
                    dict_name = item.name+'_pore_label_'+key
                    self._pore_data_amalgamate.update({dict_name : item._pore_info[key]})
        #Add geometry data
        for key in self._pore_data.keys():
            if sp.amax(self._pore_data[key]) < sp.inf:
                dict_name = 'pore'+'_'+key
                self._pore_data_amalgamate.update({dict_name : self._pore_data[key]})
        for key in self._pore_info.keys():
            if sp.amax(self._pore_info[key]) < sp.inf:
                dict_name = 'pore'+'_label_'+key
                self._pore_data_amalgamate.update({dict_name : self._pore_info[key]})
        return self._pore_data_amalgamate

    def amalgamate_throat_data(self,fluids='all'):
        r"""
        Returns a dictionary containing ALL throat data from all fluids, physics and geometry objects
        """
        self._throat_data_amalgamate = {}
        if type(fluids)!= sp.ndarray and fluids=='all':
            fluids = self._fluids
        elif type(fluids)!= sp.ndarray: 
            fluids = sp.array(fluids,ndmin=1)
        #Add fluid data
        for item in fluids:
            if type(item)==sp.str_: item =  self.find_object_by_name(item)
            for key in item._throat_data.keys():
                if sp.amax(item._throat_data[key]) < sp.inf:
                    dict_name = item.name+'_throat_'+key
                    self._throat_data_amalgamate.update({dict_name : item._throat_data[key]})
            for key in item._throat_info.keys():
                if sp.amax(item._throat_info[key]) < sp.inf:
                    dict_name = item.name+'_throat_label_'+key
                    self._throat_data_amalgamate.update({dict_name : item._throat_info[key]})
        #Add geometry data
        for key in self._throat_data.keys():
            if sp.amax(self._throat_data[key]) < sp.inf:
                dict_name = 'throat'+'_'+key
                self._throat_data_amalgamate.update({dict_name : self._throat_data[key]})
        for key in self._throat_info.keys():
            if sp.amax(self._throat_info[key]) < sp.inf:
                dict_name = 'throat'+'_label_'+key
                self._throat_data_amalgamate.update({dict_name : self._throat_info[key]})
        return self._throat_data_amalgamate

    #--------------------------------------------------------------------------
    '''Graph theory and topology related methods'''
    #--------------------------------------------------------------------------
    def create_adjacency_matrix(self,data=None,prop=None,sprsfmt='all',dropzeros=True,sym=True):
        r"""
        Generates adjacency matricies in various sparse storage formats

        Parameters
        ----------
        data : array_like (optional)
            An ndarray containing the throat values to place into the I,J locations of the IJV sparse matrix.
            If no argument is supplied then the standard adjacency matrix is assumed.
        prop : String (optional)
            The name of the property being written to the adjacency matrix.
            If no argument is supplied then the standard adjacency matrix is assumed.
        sprsfmt : String, optional
            The sparse storage format to use. If no type is specified then all are generated (coo, csr & lil)
        dropzeros : Boolean, optional
            Remove 0 elements from the values, instead of creating 0-weighted links, the default is True.
        sym : Boolean, optional
            Makes the matrix symmetric about the diagonal, the default is true.

        Returns
        -------
        adj_mat : sparse_matrix
            Returns adjacency matrix in specified format for local use.

        Notes
        -----
        This also stores the adjacency matrix in a nested dictionary called
        _adjacency_matrix.  This top level keys are the storage type, and
        under each type they keys are the property name.

        Examples
        --------
        >>> pn = OpenPNM.Network.Cubic(name='doc_test').generate(divisions=[5,5,5],lattice_spacing=[1])
        >>> vals = pn.get_throat_data(prop='numbering')
        >>> temp = pn.create_adjacency_matrix(data=vals,prop='numbering',sprsfmt='csr')
        >>> print(pn.adjacency_matrix['csr'].keys())
        dict_keys(['numbering'])

        """
        self._logger.debug('create_adjacency_matrix: Start of method')
        Np   = self.num_pores()
        Nt   = self.num_throats()

        if (data is not None) and (prop is not None):
            dataset = data
            tprop = prop
            if sp.shape(dataset)[0] != Nt:
                raise Exception('Received dataset of incorrect length')
        else:
            dataset = sp.ones(Nt)
            tprop = 'connections'

        if dropzeros:
            ind = dataset>0
        else:
            ind = sp.ones_like(dataset,dtype=bool)

        conn = self._throat_data["connections"][ind]
        row  = conn[:,0]
        col  = conn[:,1]
        data = dataset[ind]

        if sym: #Append row & col to each other, and data to itself
            row  = sp.append(row,conn[:,1])
            col  = sp.append(col,conn[:,0])
            data = sp.append(data,data)

        temp = sprs.coo_matrix((data,(row,col)),(Np,Np))
        if sprsfmt == 'coo' or sprsfmt == 'all':
            self._adjacency_matrix['coo'][tprop] = temp
        if sprsfmt == 'csr' or sprsfmt == 'all':
            self._adjacency_matrix['csr'][tprop] = temp.tocsr()
        if sprsfmt == 'lil' or sprsfmt == 'all':
            self._adjacency_matrix['lil'][tprop] = temp.tolil()
        if sprsfmt != 'all':
            return self._adjacency_matrix[sprsfmt][tprop]

    def create_incidence_matrix(self,data=None,prop=None,sprsfmt='all',dropzeros=True):
        r"""
        Creates an incidence matrix filled with supplied throat values

        Parameters
        ----------
        data : array_like (optional)
            An ndarray containing the throat values to place into the I,J locations of the IJV sparse matrix.
            If no argument is supplied then the standard adjacency matrix is assumed.
        prop : String (optional)
            The name of the property being written to the adjacency matrix.
            If no argument is supplied then the standard adjacency matrix is assumed.
        sprsfmt : String, optional
            The sparse storage format to use. If none type is given, all are generated (coo, csr & lil)
        dropzeros : Boolean, optional
            Remove 0 elements from values, instead of creating 0-weighted links, the default is True.

        Returns
        -------
        An incidence matrix (a cousin to the adjacency matrix, useful for finding throats of given a pore)

        Notes
        -----
        This also stores the incidence matrix in a nested dictionary called
        _incidence_matrix.  This top level keys are the storage type, and
        under each type they keys are the property name.

        Examples
        --------
        >>> print('nothing yet')
        nothing yet
        """
        self._logger.debug('create_incidence_matrix: Start of method')

        Nt = self.num_throats()
        Np = self.num_pores()

        if (data is not None) and (prop is not None):
            dataset = data
            tprop = prop
        else:
            dataset = sp.ones(Nt)
            tprop = 'connections'

        if dropzeros:
            ind = dataset > 0
        else:
            ind = sp.ones_like(dataset, dtype=bool)

        conn = self._throat_data['connections'][ind]
        row  = conn[:,0]
        row = sp.append(row,conn[:,1])
        col = self.get_throat_indices('all')[ind]
        col = sp.append(col,col)
        data = sp.append(dataset[ind],dataset[ind])

        temp = sprs.coo.coo_matrix((data,(row,col)),(Np,Nt))
        if sprsfmt == 'coo' or sprsfmt == 'all':
            self._incidence_matrix['coo'][tprop] = temp
        if sprsfmt == 'csr' or sprsfmt == 'all':
            self._incidence_matrix['csr'][tprop] = temp.tocsr()
        if sprsfmt == 'lil' or sprsfmt == 'all':
            self._incidence_matrix['lil'][tprop] = temp.tolil()
        if sprsfmt != 'all':
            return self._incidence_matrix[sprsfmt][tprop]
            
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
        >>> pn = OpenPNM.Network.Cubic(name='doc_test').generate(divisions=[5,5,5],lattice_spacing=[1])
        >>> pn.find_connected_pores(throats=[0,1])
        array([[0, 1],
               [0, 5]])
        >>> pn.find_connected_pores(throats=[0,1],flatten=True)
        array([0, 1, 5])
        """
        Ps = self._throat_data['connections'][throats]
        #Ps = [sp.asarray(x) for x in Ps if x]
        if flatten:
            Ps = sp.unique(sp.hstack(Ps))
        return Ps

    def find_connecting_throat(self,P1,P2):
        r"""
        Return a the throat number connecting two given pores connected

        Parameters
        ----------
        P1 , P2 : int
            The pore numbers connected by the desired throat

        Returns
        -------
        Tnum : int
            Returns throat number, or empty array if pores are not connected
            
        Examples
        --------
        >>> pn = OpenPNM.Network.Cubic(name='doc_test').generate(divisions=[5,5,5],lattice_spacing=[1])
        >>> pn.find_connecting_throat(0,1)
        array([0])
        """
        return sp.intersect1d(self.find_neighbor_throats(P1),self.find_neighbor_throats(P2))

    def find_neighbor_pores(self,pores,mode='union',flatten=True,excl_self=False):
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
        excl_self : bool, optional
            If this is True (default) then the input pores are not included
            in the returned list.  This option only applies when input pores
            are in fact neighbors to each other, otherwise they are not
            part of the returned list.  
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
        >>> pn.find_neighbor_pores(pores=[0,1]) #Find all neighbors, excluding selves (default behavior)
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
        #Count neighboring pores
        try:
            neighborPs = self._adjacency_matrix['lil']['connections'].rows[[pores]]
        except:
            self._logger.info('Creating adjacency matrix, please wait')
            self.create_adjacency_matrix()
            neighborPs = self._adjacency_matrix['lil']['connections'].rows[[pores]]
        if flatten:
            #All the empty lists must be removed to maintain data type after hstack (numpy bug?)
            neighborPs = [sp.asarray(x) for x in neighborPs if x]
            neighborPs = sp.hstack(neighborPs)
            #neighborPs = sp.concatenate((neighborPs,pores))
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
        >>> pn = OpenPNM.Network.Cubic(name='doc_test').generate(divisions=[5,5,5],lattice_spacing=[1])
        >>> pn.find_neighbor_throats(pores=[0,1])
        array([0, 1, 2, 3, 4, 5])
        >>> pn.find_neighbor_throats(pores=[0,1],flatten=False)
        array([array([0, 1, 2]), array([0, 3, 4, 5])], dtype=object)
        """
        #Test for existance of incidence matrix
        try:
            neighborTs = self._incidence_matrix['lil']['connections'].rows[[pores]]
        except:
            self._logger.info('Creating incidence matrix, please wait')
            self.create_incidence_matrix()
            neighborTs = self._incidence_matrix['lil']['connections'].rows[[pores]]
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

    def show_boundaries(self,Ptype=[]):
        r'''
        Documentation for this method is being updated, we are sorry for the inconvenience.
        '''
        pass
        #        from mpl_toolkits.mplot3d import Axes3D
        #        #Parse Ptype input argument
        #        if Ptype == [] or Ptype == 'all':
        #            Ptype = self.get_type_definitions()[1]
        #        elif type(Ptype[0]) == str:
        #            Ptype = self.get_type_definitions(Ptype)
        #        fig = plt.figure()
        #        ax = fig.add_subplot(111, projection='3d')
        #        for i in Ptype:
        #
        #            xs = self._pore_data['coords'][:,0]
        #            ys = self._pore_data['coords'][:,1]
        #            zs = self._pore_data['coords'][:,2]
        #            ax.scatter(xs, ys, zs, zdir='z', s=20, c='b')
        #        plt.show()
        
    def clone_pores(self,pores,mode='parent',apply_label=['clone']):
        r'''
        mode options should be 'parent', 'siblings'
        '''
        if sp.shape(self.props(pores='all'))[0] > 1:
            raise Exception('Cannot clone an active network')
        apply_label = list(apply_label)
        #Clone pores
        Np = self.num_pores()
        parents = sp.array(pores,ndmin=1)
        pcurrent = self.get_pore_data(prop='coords')
        pclone = pcurrent[pores,:]
        pnew = sp.concatenate((pcurrent,pclone),axis=0)
        Npnew = sp.shape(pnew)[0]
        clones = sp.arange(Np,Npnew)
        #Increase size of 'all' to accomodate new pores
        self.set_pore_info(label='all', locations=sp.ones((Npnew,),dtype=bool))
        #Insert cloned pore coordinates into network
        self.set_pore_data(prop='coords',data=pnew)
        #Apply specified labels to cloned pores
        for item in apply_label:
            self.set_pore_info(label=item,locations=clones)

        #Add connections between parents and clones
        tcurrent = self.get_throat_data(prop='connections')
        tclone = sp.vstack((parents,clones)).T
        tnew = sp.concatenate((tcurrent,tclone),axis=0)
        Ntnew = sp.shape(tnew)[0]
        #Increase size of 'all' to accomodate new throats
        self.set_throat_info(label='all', locations=sp.ones((Ntnew,),dtype=bool))
        #Insert new throats into network
        self.set_throat_data(prop='connections',data=tnew)
        
        # Any existing adjacency and incidence matrices will be invalid
        self.reset_graphs()
        
    def __str__(self):
        r"""
        Print some basic properties
        """
        self._logger.debug("Method: print_overview")

        str_overview = (
        '==================================================\n'
        'Overview of network properties\n'
        '--------------------------------------------------\n'
        'Basic properties of the network\n'
        '- Number of pores:   {num_pores}\n'
        '- Number of throats: {num_throats}\n'
        ).format(num_pores=self.num_pores(),
                   num_throats=self.num_throats())

        str_pore = "\nPore properties:"
        for key, value in self._pore_data.items():
            str_pore += "\n\t{0:20}{1.dtype:20}{1.shape:20}".format(key,value)

        str_throat = "\nThroat properties:"
        for key, value in self._throat_data.items():
            str_throat += "\n\t{0:20}{1.dtype:20}{1.shape:20}".format(key,value)


        return str_overview+str_pore+str_throat

    def regenerate_fluids(self):
        r'''
        '''
        fluids = []
        for item in self._fluids.keys():
            fluids.append(self._fluids[item])
        for item in fluids:
            self._logger.info('Regenerating properties for '+item.name)
            for prop in item._prop_list:
                item.regenerate(prop_list=prop)

    def regenerate_physics(self):
        r'''
        '''
        physics = []
        for item1 in self._fluids.keys():
            for item2 in self._fluids[item1]._physics.keys():
                physics.append(self._fluids[item1]._physics[item2])
        for item in physics:
            self._logger.info('Regenerating pore properties for '+item.name)
            for prop in item._prop_list:
                item.regenerate(prop_list=prop)
                
    def regenerate_geometries(self):
        r'''
        '''
        geometry = []
        for item in self._geometries.keys():
            geometry.append(self._geometries[item])
        #Regenerate pores first
        for item in geometry:
            self._logger.info('Regenerating pore properties for '+item.name)
            for prop in item._prop_list:
                if prop.split('_')[0] == 'pore':
                    item.regenerate(prop_list=prop)
        #Regenerate throats second
        for item in geometry:
            self._logger.info('Regenerating throat properties for '+item.name)
            for prop in item._prop_list:
                if prop.split('_')[0] == 'throat':
                    item.regenerate(prop_list=prop)
            
    def add_geometry(self,name=None,subclass='GenericGeometry',**kwargs):
        r'''
        '''
        temp = OpenPNM.Geometry.__getattribute__(subclass)
        return temp(network=self,name=name,**kwargs)

    def add_fluid(self,name=None,subclass='GenericFluid',**kwargs):
        r'''
        '''
        temp = OpenPNM.Fluids.__getattribute__(subclass)
        return temp(network=self,name=name,**kwargs)
        
    def add_physics(self,fluid,geometry,name=None,subclass='GenericPhysics',**kwargs):
        r'''
        '''
        temp = OpenPNM.Physics.__getattribute__(subclass)
        return temp(network=self,name=name,fluid=fluid,geometry=geometry,**kwargs)

    def update_network(self):
        r'''
        Documentation for this method is being updated, we are sorry for the inconvenience.
        '''
        self.create_adjacency_matrix()
        self.create_incidence_matrix()
        
    def reset_graphs(self):
        r'''
        Clears the adjacency and incidence matrices
        '''
        #Re-initialize adjacency and incidence matrix dictionaries
        self._logger.debug('Resetting adjacency and incidence matrices')
        self._adjacency_matrix['coo'] = {}
        self._adjacency_matrix['csr'] = {}
        self._adjacency_matrix['lil'] = {}
        self._incidence_matrix['coo'] = {}
        self._incidence_matrix['csr'] = {}
        self._incidence_matrix['lil'] = {}        

    def trim(self, pores=[], throats=[]):
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
        self.get_pore_data(prop='diameter') < 0.5 # [True, False, ...]
        And send it over as an argument. 
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.trim(pores=[35])
        '''
        pores = np.ravel(pores)
        throats = np.ravel(throats)
        
        if sp.shape(pores)[0]>0:
            Pdrop = sp.zeros((self.num_pores(),),dtype=bool)
            Pdrop[pores] = True
            Pkeep = ~Pdrop
            Tdrop = sp.zeros((self.num_throats(),),dtype=bool)
            Ts = self.find_neighbor_throats(pores)
            Tdrop[Ts] = 1
            Tkeep = ~Tdrop
        if sp.shape(throats)[0]>0:
            Tdrop = sp.zeros((self.num_throats(),),dtype=bool)
            Tdrop[throats] = 1
            Tkeep = ~Tdrop
            Pkeep = self.get_pore_indices(labels='all')
            Pkeep = self.to_mask(pores=Pkeep)
        
        #Remap throat connections
        Pnew = sp.arange(0,sum(Pkeep),dtype=int)
        Pmap = sp.ones((self.num_pores(),),dtype=int)*-1
        Pmap[Pkeep] = Pnew
        tpore1 = self.get_throat_data(prop='connections')[:,0]
        tpore2 = self.get_throat_data(prop='connections')[:,1]
        temp1 = Pmap[tpore1[Tkeep]]
        temp2 = Pmap[tpore2[Tkeep]]
        
        # Insert new indices into network
        # Write connections specifically
        self._throat_data['connections'] = sp.vstack((temp1,temp2)).T
        # Over-write remaining data
        for item in self._pore_data.keys():
            self._pore_data[item] = self._pore_data[item][Pkeep]
        for item in self._pore_info.keys():
            self._pore_info[item] = self._pore_info[item][Pkeep]            
        for item in self._throat_data.keys():
            if item != 'connections':
                self._throat_data[item] = self._throat_data[item][Tkeep]  
        for item in self._throat_info.keys():
            self._throat_info[item] = self._throat_info[item][Tkeep]
        
        #Reset network
        self.reset_graphs()
        
        #Check for individual isolated pores
        Ps = sp.sum(self.num_neighbors(self.get_pore_indices())==0)
        if Ps > 0:
            self._logger.warning(str(Ps)+' pores no longer have neighbors')
        
        #Check for clusters of isolated pores
        Cs = self.find_clusters(self.throats('all'))
        if sp.shape(sp.unique(Cs))[0] > (Ps+1):
            self._logger.warning('Isolated clusters exist in the network')
        
    def find_clusters(self,mask):
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
        self.create_adjacency_matrix(prop='temp', data=temp, sprsfmt='csr', dropzeros=True)
        clusters = sprs.csgraph.connected_components(self._adjacency_matrix['csr']['temp'])[1]
        del self._adjacency_matrix['csr']['temp']
        return clusters

if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)
    
    
    
    
    
    
    
    
    
    
