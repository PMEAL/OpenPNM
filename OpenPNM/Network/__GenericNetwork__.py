"""
module __GenericNetwork__: Abstract class to construct pore networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Network.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, parent_dir)
import OpenPNM

import scipy as sp
import scipy.sparse as sprs
import pprint

class GenericNetwork(OpenPNM.Base.Tools):
    r"""
    GenericNetwork - Base class to construct pore networks

    This class contains the interface definition for the construction of networks

    Parameters
    ----------
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
        self._logger.debug("Construct Network container")
        self.name = name
        #Initialize adjacency and incidence matrix dictionaries
        self.adjacency_matrix = {}
        self.incidence_matrix = {}
        self.adjacency_matrix['coo'] = {}
        self.adjacency_matrix['csr'] = {}
        self.adjacency_matrix['lil'] = {}
        self.incidence_matrix['coo'] = {}
        self.incidence_matrix['csr'] = {}
        self.incidence_matrix['lil'] = {}
        self._logger.info("Construction of Network container complete")
        
    def generate(self, **params):
        r"""
        Generate the network
        """
        print('this method is not implemented')
 
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

        Notes
        -----
        This uses an unweighted average, without attempting to account for distances or sizes of pores and throats.

        """
        if sp.size(Tvals)==1:
            Pvals = Tvals
        elif sp.size(Tvals) != self.get_num_throats():
            raise Exception('The list of throat information received was the wrong length')
        else:
            Pvals = sp.zeros((self.get_num_pores()))
            #Only interpolate conditions for internal pores, type=0
            Pnums = sp.r_[0:self.get_num_pores(Ptype=[0])]
            nTs = self.get_neighbor_throats(Pnums,flatten=False)
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

        Notes
        -----
        This uses an unweighted average, without attempting to account for distances or sizes of pores and throats.

        """
        if sp.size(Pvals)==1:
            Tvals = Pvals
        elif sp.size(Pvals) != self.get_num_pores():
            raise Exception('The list of pore information received was the wrong length')
        else:
            Tvals = sp.zeros((self.get_num_throats()))
            #Interpolate values for all throats, including those leading to boundary pores
            Tnums = sp.r_[0:self.get_num_throats()]
            nPs = self.get_connected_pores(Tnums,flatten=False)
            for i in sp.r_[0:sp.shape(nPs)[0]]:
                Tvals[i] = sp.mean(Pvals[nPs[i]])
        return Tvals

    def amalgamate_pore_data(self):
        r"""
        Returns a dictionary containing ALL pore data from all fluids, physics and geometry objects
        """
        self._pore_data = {}
        #Add fluid conditions
        for item in self._fluids:
            for key in item.pore_data.keys():
                dict_name = item.name+'_'+key
                self._pore_data.update({dict_name : item.pore_data[key]})
        #Add physics conditions (does nothing now since they're stored in fluids)
        for item in self._physics:
            for key in item.pore_data.keys():
                dict_name = item.name+'_'+key
                self._pore_data.update({dict_name : item.pore_data[key]})
        #Add geometry data
        for key in self._pore_data.keys():
            dict_name = 'pore'+'_'+key
            self._pore_data.update({dict_name : self._pore_data[key]})
        pprint.pprint(self._pore_data)
        return self._pore_data

    def amalgamate_throat_data(self):
        r"""
        Returns a dictionary containing ALL throat data from all fluids, physics and geometry objects
        """
        self._throat_data = {}
        #Add fluid conditions
        for item in self._fluids:
            for key in item.throat_data.keys():
                dict_name = item.name+'_'+key
                self._throat_data.update({dict_name : item.throat_data[key]})
        #Add physics conditions (does nothing now since they're stored in fluids)
        for item in self._physics:
            for key in item.throat_data.keys():
                dict_name = item.name+'_'+key
                self._throat_data.update({dict_name : item.throat_conditions[key]})
        #Add geometry data
        for key in self._throat_data.keys():
            dict_name = 'throat'+'_'+key
            self._throat_data.update({dict_name : self._throat_data[key]})
        pprint.pprint(self._throat_data)
        return self._throat_data

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
        Np   = self.get_num_pores()
        Nt   = self.get_num_throats()

        if (data!=None) and (prop!=None):
            dataset = data
            tprop = prop
            if sp.shape(dataset)[0]!=Nt:
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

        Nt = self.get_num_throats()
        Np = self.get_num_pores()

        if (data!=None) and (prop!=None):
            dataset = data
            tprop = prop
        else:
            dataset = sp.ones(Nt)
            tprop = 'connections'

        if dropzeros:
            ind = dataset>0
        else:
            ind = sp.ones_like(dataset,dtype=bool)

        conn = self._throat_data['connections'][ind]
        row  = conn[:,0]
        row = sp.append(row,conn[:,1])
        col = self._throat_data['numbering'][ind]
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

    def get_connected_pores(self,tnums=[],flatten=False):
        r"""
        Return a list of pores connected to a list of throats

        Parameters
        ----------
        Tnums : array_like
            List of throats ID numbers
        flatten : boolean, optional
            If flatten is True (default) a 1D array of unique pore ID numbers
            is returned. If flatten is False the returned array contains
            arrays of neighboring pores for each input throat, in the order
            they were sent.

        Returns
        -------
        Ps : 1D array (if flatten is True) or ndarray of arrays (if flatten is
            False)

        Examples
        --------
        >>> pn = OpenPNM.Network.Cubic(name='doc_test').generate(divisions=[5,5,5],lattice_spacing=[1])
        >>> pn.get_connected_pores(tnums=[0,1])
        array([[0, 1],
               [0, 5]])
        >>> pn.get_connected_pores(tnums=[0,1],flatten=True)
        array([0, 1, 5])
        """
        Ps = self._throat_data['connections'][tnums]
        #Ps = [sp.asarray(x) for x in Ps if x]
        if flatten:
            Ps = sp.unique(sp.hstack(Ps))
        return Ps

    def get_connecting_throat(self,P1,P2):
        r"""
        Return a the throat number connecting two given pores connected

        Parameters
        ----------
        Pnum1 , Pnum2 : int
            The pore numbers connected by the desired throat

        Returns
        -------
        Tnum : int
            Returns throat ID number, or empty array if pores are not connected
            
        Examples
        --------
        >>> pn = OpenPNM.Network.Cubic(name='doc_test').generate(divisions=[5,5,5],lattice_spacing=[1])
        >>> pn.get_connecting_throat(0,1)
        array([0])
        """
        return sp.intersect1d(self.get_neighbor_throats(P1),self.get_neighbor_throats(P2))

    def get_neighbor_pores(self,pnums,subdomain=['all'],flatten=True,mode=''):
        r"""
        Returns a list of pores neighboring the given pore(s)

        Parameters
        ----------
        pnums : array_like
            ID numbers of pores whose neighbors are sought.
        subdomain : list of strings
            Type of pores to be returned
        flatten : boolean, optional
            If flatten is True (default) a 1D array of unique pore ID numbers
            is returned with the input pores (Pnum) removed. If flatten is
            False the returned array contains arrays of neighboring pores for
            each input pore, in the order they were sent.
        mode : string, optional
            This controls the contents of what indices are returned.  The mode can only be specified to a 'flattened' list.
            Options are 'union', 'intersection', and 'not_intersection'

        Returns
        -------
        neighborPs : 1D array (if flatten is True) or ndarray of ndarrays (if flatten if False)

        Examples
        --------
        >>> pn = OpenPNM.Network.Cubic(name='doc_test').generate(divisions=[5,5,5],lattice_spacing=[1])
        >>> pn.get_neighbor_pores(pnums=[0,1])
        array([ 2,  5,  6, 25, 26])
        >>> pn.get_neighbor_pores(pnums=[0,1],flatten=False)
        array([array([ 1,  5, 25]), array([ 0,  2,  6, 26])], dtype=object)
        """
        #Convert string to list, if necessary
        if type(subdomain) == str: subdomain = [subdomain]
        #Count neighboring pores
        try:
            neighborPs = self.adjacency_matrix['lil']['connections'].rows[[pnums]]
        except:
            self._logger.info('Creating adjacency matrix, please wait')
            self.create_adjacency_matrix()
            neighborPs = self.adjacency_matrix['lil']['connections'].rows[[pnums]]
        if flatten:
            #All the empty lists must be removed to maintain data type after hstack (numpy bug?)
            neighborPs = [sp.asarray(x) for x in neighborPs if x]
            neighborPs = sp.hstack(neighborPs)
            #neighborPs = sp.concatenate((neighborPs,pnums))
            #Remove references to input pores and duplicates
            if mode == 'not_intersection':
                neighborPs = sp.unique(sp.where(sp.bincount(neighborPs)==1)[0])
            elif mode == 'union':
                neighborPs = sp.unique(neighborPs)
            elif mode == 'intersection':
                neighborPs = sp.unique(sp.where(sp.bincount(neighborPs)>1)[0])
            elif mode == '':
                neighborPs = sp.unique(neighborPs[~sp.in1d(neighborPs,pnums)])
            #Remove pores of the wrong type
            mask = self.get_pore_indices(subdomain=subdomain,indices=False)
            neighborPs = neighborPs[mask[neighborPs]]
        else:
            mask = self.get_pore_indices(subdomain=subdomain,indices=False)
            for i in range(0,sp.size(pnums)):
                neighborPs[i] = sp.array(neighborPs[i])[mask[neighborPs[i]]]
        return sp.array(neighborPs,ndmin=1)

    def get_neighbor_throats(self,pnums,subdomain=['all'],flatten=True,mode=''):
        r"""
        Returns a list of throats neighboring the given pore(s)

        Parameters
        ----------
        pnums : array_like
            ID numbers of pores whose neighbors are sought
        flatten : boolean, optional
            If flatten is True (default) a 1D array of unique throat ID numbers
            is returned. If flatten is False the returned array contains arrays
            of neighboring throat ID numbers for each input pore, in the order
            they were sent.
        mode : string, optional
            This controls the contents of what indices are returned.  The mode can only be specified to a 'flattened' list.
            Options are 'union', 'intersection', and 'not_intersection'

        Returns
        -------
        neighborTs : 1D array (if flatten is True) or ndarray of arrays (if
            flatten if False)

        Examples
        --------
        >>> pn = OpenPNM.Network.Cubic(name='doc_test').generate(divisions=[5,5,5],lattice_spacing=[1])
        >>> pn.get_neighbor_throats(pnums=[0,1])
        array([0, 1, 2, 3, 4, 5])
        >>> pn.get_neighbor_throats(pnums=[0,1],flatten=False)
        array([array([0, 1, 2]), array([0, 3, 4, 5])], dtype=object)
        """
        #Convert string to list, if necessary
        if type(subdomain) == str: subdomain = [subdomain]
        #Test for existance of incidence matrix
        try:
            neighborTs = self.incidence_matrix['lil']['connections'].rows[[pnums]]
        except:
            self._logger.info('Creating incidence matrix, please wait')
            self.create_incidence_matrix()
            neighborTs = self.incidence_matrix['lil']['connections'].rows[[pnums]]
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
            elif mode == '':
                neighborTs = sp.unique(neighborTs)
            #Remove throats of the wrong type            
            mask = self.get_throat_indices(subdomain=subdomain,indices=False)
            neighborTs = neighborTs[mask[neighborTs]]
        else:
            mask = self.get_throat_indices(subdomain=subdomain,indices=False)
            for i in range(0,sp.size(pnums)):
                neighborTs[i] = sp.array(neighborTs[i])[mask[neighborTs[i]]]
        return sp.array(neighborTs,ndmin=1)

    def get_num_neighbors(self,pnums,subdomain=['all'],flatten=True):
        r"""
        Returns an ndarray containing the number of pores for each element in Pnums

        Parameters
        ----------
        pnums : array_like
            ID numbers of pores whose neighbors are sought

        Returns
        -------
        num_neighbors : 1D array with number of neighbors in each element

        Examples
        --------
        >>> pn = OpenPNM.Network.Cubic(name='doc_test').generate(divisions=[5,5,5],lattice_spacing=[1])
        >>> Pnum = [0,1]
        >>> pn.get_num_neighbors(Pnum,flatten=False)
        array([3, 4], dtype=int8)
        >>> pn.get_num_neighbors(Pnum)
        7
        """
        #Convert string to list, if necessary
        if type(subdomain) == str: subdomain = [subdomain]

        #Count number of neighbors
        neighborPs = self.get_neighbor_pores(pnums,subdomain=subdomain,flatten=False)
        num = sp.zeros(sp.shape(neighborPs),dtype=sp.int8)
        for i in range(0,sp.shape(num)[0]):
            num[i] = sp.size(neighborPs[i])
        if flatten:
            num = sp.sum(num)
        return num

    def check_basic(self):
        r"""
        Check the network for general health
        TODO: implement
        """
        self._logger.debug("Method: check for general healts")

    def print_overview(self):
        r"""
        Print some basic properties
        """
        pprint.pprint(self._pore_data)
        pprint.pprint(self._throat_data)

    def show_boundaries(self,Ptype=[]):
        r'''
        '''
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
        ).format(num_pores=self.get_num_pores(),
                   num_throats=self.get_num_throats())

        str_pore = "\nPore properties:"
        for key, value in self._pore_data.iteritems():
            str_pore += "\n\t{0:20}{1.dtype:20}{1.shape:20}".format(key,value)

        str_throat = "\nThroat properties:"
        for key, value in self._throat_data.iteritems():
            str_throat += "\n\t{0:20}{1.dtype:20}{1.shape:20}".format(key,value)

#        str_pore_cond = "\nPore conditions:"
#        for key, value in self.pore_conditions.iteritems():
#            str_pore_cond += "\n\t{0:20}{1.dtype:20}{1.shape:20}".format(key,sp.array(value))
#
#        str_throat_cond = "\nThroat conditions:"
#        for key, value in self.throat_conditions.iteritems():
#            str_throat_cond += "\n\t{0:20}{1.dtype:20}{1.shape:20}".format(key,sp.array(value))

        return str_overview+str_pore+str_throat

    def update_network(self):
        r"""
        """
        self.create_adjacency_matrix()
        self.create_incidence_matrix()

if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)
