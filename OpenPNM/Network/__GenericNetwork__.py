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
import scipy as sp
import scipy.sparse as sprs
import pprint

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

    def __init__(self, name,**kwargs):
        r"""
        Initialize
        """
        super(GenericNetwork,self).__init__(**kwargs)
        self._logger.info("Construct Network")
        self.name = name
        #Initialize fluid, physics, and geometry tracking lists
        self._fluids = {}
        self._geometries ={}
        self._physics = {}
        #Initialize adjacency and incidence matrix dictionaries
        self.adjacency_matrix = {}
        self.incidence_matrix = {}
        self.adjacency_matrix['coo'] = {}
        self.adjacency_matrix['csr'] = {}
        self.adjacency_matrix['lil'] = {}
        self.incidence_matrix['coo'] = {}
        self.incidence_matrix['csr'] = {}
        self.incidence_matrix['lil'] = {}
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
            Pnums = sp.r_[0:self.num_pores(Ptype=[0])]
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
        adj_mat : sparse_matrix, optional
            Returns adjacency matrix in specified format for local use.

        Notes
        -----
        This 'can' return the specified sparse matrix, but will always write the generated matrix to the network object

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

        #Append row & col to each other, and data to itself
        if sym:
            row  = sp.append(row,conn[:,1])
            col  = sp.append(col,conn[:,0])
            data = sp.append(data,data)

        temp = sprs.coo_matrix((data,(row,col)),(Np,Np))
        if sprsfmt == 'coo' or sprsfmt == 'all':
            self.adjacency_matrix['coo'][tprop] = temp
        if sprsfmt == 'csr' or sprsfmt == 'all':
            self.adjacency_matrix['csr'][tprop] = temp.tocsr()
        if sprsfmt == 'lil' or sprsfmt == 'all':
            self.adjacency_matrix['lil'][tprop] = temp.tolil()
        if sprsfmt != 'all':
            return self.adjacency_matrix[sprsfmt][tprop]

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
            self.incidence_matrix['coo'][tprop] = temp
        if sprsfmt == 'csr' or sprsfmt == 'all':
            self.incidence_matrix['csr'][tprop] = temp.tocsr()
        if sprsfmt == 'lil' or sprsfmt == 'all':
            self.incidence_matrix['lil'][tprop] = temp.tolil()
        if sprsfmt != 'all':
            return self.incidence_matrix[sprsfmt][tprop]
        
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
        
    def clone_pores(self,pnums,mode='parent',apply_label=['clone']):
        r'''
        mode options should be 'parent', 'siblings'
        '''
        if type(pnums) == str: 
            pnums = self.get_pore_indices(labels=[pnums])
        if self._geometries != {} or self._fluids != {}:
            raise Exception('Cannot clone an active network')
        apply_label = list(apply_label)
        #Clone pores
        Np = self.num_pores()
        parents = sp.array(pnums,ndmin=1)
        pcurrent = self.get_pore_data(prop='coords')
        pclone = pcurrent[pnums,:]
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
        self._reset_network()
        
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

    def regenerate_fluid(self,fluids='all',prop_list=''):
        
        if type(fluids)!= sp.ndarray and fluids=='all':
            fluids = self._fluids
        elif type(fluids)!= sp.ndarray: 
            fluids = sp.array(fluids,ndmin=1)
        for item in fluids:
            try: item = self.find_object_by_name(item) 
            except: pass #Accept object
            item.regenerate(prop_list=prop_list)
      

    def regenerate_physics(self,physics='all',prop_list=''):
        
        if type(physics)!= sp.ndarray and physics=='all':
            physics = self._physics
        elif type(physics)!= sp.ndarray: 
            physics = sp.array(physics,ndmin=1)
        for item in physics:
            try: item = self.find_object_by_name(item) 
            except: pass #Accept object
            item.regenerate(prop_list=prop_list) 
                
    def regenerate_geometry(self,geometry='all',prop_list=''):
        
        if type(geometry)!= sp.ndarray and geometry=='all':
            geometry = self._geometry
        elif type(geometry)!= sp.ndarray: 
            geometry = sp.array(geometry,ndmin=1)
        for item in geometry:
            try: item = self.find_object_by_name(item) 
            except: pass #Accept object
            item.regenerate(prop_list=prop_list) 
            
    def add_geometry(self,name,subclass='GenericGeometry',**kwargs):
        r'''
        '''
        temp = OpenPNM.Geometry.__getattribute__(subclass)
        return temp(network=self,name=name,**kwargs)

    def add_fluid(self,name,subclass='GenericFluid',**kwargs):
        r'''
        '''
        temp = OpenPNM.Fluids.__getattribute__(subclass)
        return temp(network=self,name=name,**kwargs)
        
    def add_physics(self,name,fluid,geometry,subclass='GenericPhysics',**kwargs):
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
        
    def _reset_network(self):
        r'''
        '''
        #Initialize adjacency and incidence matrix dictionaries
        self.adjacency_matrix = {}
        self.incidence_matrix = {}
        self.adjacency_matrix['coo'] = {}
        self.adjacency_matrix['csr'] = {}
        self.adjacency_matrix['lil'] = {}
        self.incidence_matrix['coo'] = {}
        self.incidence_matrix['csr'] = {}
        self.incidence_matrix['lil'] = {}
        
    def pores(self,label='all'):
        r'''
        Returns a list of pore indices
        
        Parameters
        ----------
        label : string
            Accepts a *single* label specifying the pore label of interest
            
        Returns
        -------
        indices : array_like
            A (1,) array of pore indices
        
        See Also
        --------
        get_pore_indices
        
        Notes
        -----
        This is a helper function to simplify the process of retrieving pores
        
        '''
        if type(label) == str:
            return self.get_pore_indices(labels=label)
        else:
            raise Exception('The pores() method only accepts a single label, for more complex queries use get_pore_indicies')
        
        
    def throats(self,label='all'):
        r'''
        Returns a list of throat indices
        
        Parameters
        ----------
        label : string
            Accepts a *single* label specifying the throat label of interest
            
        Returns
        -------
        indices : array_like
            A (1,) array of throat indices
        
        See Also
        --------
        get_throat_indices
        
        Notes
        -----
        This is a helper function to simplify the process of retrieving throats
        
        '''
        if type(label) == str:
            return self.get_throat_indices(labels=label)
        else:
            raise Exception('The throats() method only accepts a single label, for more complex queries use get_throat_indicies')
        
        
        
        

if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)
    
    
    
    
    
    
    
    
    
    
