# -*- coding: utf-8 -*-
"""
===============================================================================
GenericNetwork: Abstract class to construct pore networks
===============================================================================

"""
import scipy as sp
import scipy.sparse as sprs
import OpenPNM.Utilities.misc as misc
from OpenPNM.Base import Core, Controller, Tools, logging
logger = logging.getLogger(__name__)
ctrl = Controller()

class GenericNetwork(Core):
    r"""
    GenericNetwork - Base class to construct pore networks

    Parameters
    ----------
    name : string
        Unique name for Network object

    """
    def __init__(self,**kwargs):
        r"""
        Initialize Network
        """
        super(GenericNetwork,self).__init__(**kwargs)
        logger.name = self.name

        #Initialize adjacency and incidence matrix dictionaries
        self._incidence_matrix = {}
        self._adjacency_matrix = {}

    def __setitem__(self,prop,value):
        for geom in self._geometries:
            if (prop in geom.keys()) and ('all' not in prop.split('.')):
                logger.error(prop+' is already defined in at least one associated Geometry object')
                return
        super(GenericNetwork,self).__setitem__(prop,value)

    def __getitem__(self,key):
        if key.split('.')[-1] == self.name:
            element = key.split('.')[0]
            return self[element+'.all']
        if key not in self.keys():
            logger.debug(key+' not on Network, constructing data from Geometries')
            return self._interleave_data(key,self.geometries())
        else:
            return super(GenericNetwork,self).__getitem__(key)
            
    def _set_net(self,network):
        pass
    def _get_net(self):
        return self
    _net = property(fset=_set_net,fget=_get_net)
    
    #--------------------------------------------------------------------------
    '''Graph Theory and Network Query Methods'''
    #--------------------------------------------------------------------------
    def create_adjacency_matrix(self,data=None,sprsfmt='coo',dropzeros=True,sym=True):
        r"""
        Generates a weighted adjacency matrix in the desired sparse format

        Parameters
        ----------
        data : array_like, optional
            An array containing the throat values to enter into the matrix (in
            graph theory these are known as the 'weights').  If omitted, ones
            are used to create a standard adjacency matrix representing
            connectivity only.

        sprsfmt : string, optional
            The sparse storage format to return.  Options are:

            * 'coo' : (default) This is the native format of OpenPNM data

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
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> vals = sp.rand(pn.num_throats(),) < 0.5
        >>> temp = pn.create_adjacency_matrix(data=vals,sprsfmt='csr')

        """
        logger.debug('create_adjacency_matrix: Start of method')
        Np   = self.num_pores()
        Nt   = self.num_throats()

        #Check if provided data is valid
        if data is None:
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
        logger.debug('create_adjacency_matrix: End of method')
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
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> vals = sp.rand(pn.num_throats(),) < 0.5
        >>> temp = pn.create_incidence_matrix(data=vals,sprsfmt='csr')
        """
        logger.debug('create_incidence_matrix: Start of method')

        Nt = self.num_throats()
        Np = self.num_pores()

        #Check if provided data is valid
        if data is None:
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
        col = self.throats('all')[ind]
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
        logger.debug('create_incidence_matrix: End of method')
        return temp
         
    def find_connected_pores(self,throats=[],flatten=False):
        r"""
        Return a list of pores connected to the given list of throats

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
        >>> import OpenPNM
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
        P1 , P2 : array_like
            The pore numbers whose throats are sought.  These can be vectors
            of pore numbers, but must be the same length

        Returns
        -------
        Tnum : list of list of int
            Returns throat number(s), or empty array if pores are not connected

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.find_connecting_throat([0,1,2],[2,2,2])
        [[], [3], []]

        TODO: This now works on 'vector' inputs, but is not actually vectorized
        in the Numpy sense, so could be slow with large P1,P2 inputs
        """
        Ts1 = self.find_neighbor_throats(P1,flatten=False)
        Ts2 = self.find_neighbor_throats(P2,flatten=False)
        Ts = []
        if sp.shape(P1) == ():
            P1 = [P1]
            P2 = [P2]
        for row in range(0,len(P1)):
            if P1[row] == P2[row]:
                throat = []
            else:
                throat = sp.intersect1d(Ts1[row],Ts2[row]).tolist()
            Ts.insert(0,throat)
        Ts.reverse()
        return Ts

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
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.find_neighbor_pores(pores=[0,2])
        array([ 1,  3,  5,  7, 25, 27])
        >>> pn.find_neighbor_pores(pores=[0,1]) #Find all neighbors, excluding selves
        array([ 2,  5,  6, 25, 26])
        >>> pn.find_neighbor_pores(pores=[0,1],mode='union',excl_self=False) #Find all neighbors, including selves
        array([ 0,  1,  2,  5,  6, 25, 26])
        >>> pn.find_neighbor_pores(pores=[0,2],flatten=False)
        array([array([ 1,  5, 25]), array([ 1,  3,  7, 27])], dtype=object)
        >>> pn.find_neighbor_pores(pores=[0,2],mode='intersection') #Find only common neighbors
        array([1])
        >>> pn.find_neighbor_pores(pores=[0,2],mode='not_intersection') #Exclude common neighbors
        array([ 3,  5,  7, 25, 27])
        """
        pores = sp.array(pores,ndmin=1)
        try:
            neighborPs = self._adjacency_matrix['lil'].rows[[pores]]
        except:
            temp = self.create_adjacency_matrix(sprsfmt='lil')
            self._adjacency_matrix['lil'] = temp
            neighborPs = self._adjacency_matrix['lil'].rows[[pores]]
        if [sp.asarray(x) for x in neighborPs if x] == []:
            return sp.array([],ndmin=1)
        if flatten:
            #All the empty lists must be removed to maintain data type after hstack (numpy bug?)
            neighborPs = [sp.asarray(x) for x in neighborPs if x]
            neighborPs = sp.hstack(neighborPs)
            neighborPs = sp.concatenate((neighborPs,pores))
            #Remove references to input pores and duplicates
            if mode == 'not_intersection':
                neighborPs = sp.array(sp.unique(sp.where(sp.bincount(neighborPs)==1)[0]),dtype=int)
            elif mode == 'union':
                neighborPs = sp.array(sp.unique(neighborPs),int)
            elif mode == 'intersection':
                neighborPs = sp.array(sp.unique(sp.where(sp.bincount(neighborPs)>1)[0]),dtype=int)
            if excl_self:
                neighborPs = neighborPs[~sp.in1d(neighborPs,pores)]
        else:
            for i in range(0,sp.size(pores)):
                neighborPs[i] = sp.array(neighborPs[i],dtype=int)
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
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.find_neighbor_throats(pores=[0,1])
        array([0, 1, 2, 3, 4, 5])
        >>> pn.find_neighbor_throats(pores=[0,1],flatten=False)
        array([array([0, 1, 2]), array([0, 3, 4, 5])], dtype=object)
        """
        #Test for existence of incidence matrix
        try:
            neighborTs = self._incidence_matrix['lil'].rows[[pores]]
        except:
            temp = self.create_incidence_matrix(sprsfmt='lil')
            self._incidence_matrix['lil'] = temp
            neighborTs = self._incidence_matrix['lil'].rows[[pores]]
        if [sp.asarray(x) for x in neighborTs if x] == []:
            return sp.array([],ndmin=1)
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
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.num_neighbors(pores=[0,1],flatten=False)
        array([3, 4])
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
            num = sp.zeros(sp.shape(neighborPs),dtype=int)
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
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn['pore.domain1'] = False
        >>> pn['pore.domain2'] = False
        >>> pn['pore.domain1'][[0,1,2]] = True
        >>> pn['pore.domain2'][[5,6,7]] = True
        >>> pn.find_interface_throats(labels=['domain1','domain2'])
        array([1, 4, 7])
        '''
        Tind = sp.array([],ndmin=1)
        if sp.shape(labels)[0] != 2:
            logger.error('Exactly two labels must be given')
            pass
        else:
            P1 = self.pores(labels=labels[0])
            P2 = self.pores(labels=labels[1])
            #Check if labels overlap
            if sp.sum(sp.in1d(P1,P2)) > 0:
                logger.error('Some labels overlap, iterface cannot be found')
                pass
            else:
                T1 = self.find_neighbor_throats(P1)
                T2 = self.find_neighbor_throats(P2)
                Tmask = sp.in1d(T1,T2)
                Tind = T1[Tmask]
        return Tind

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
        clusters = sprs.csgraph.connected_components(csgraph=temp,directed=False)[1]
        return clusters

    #--------------------------------------------------------------------------
    '''Network Manipulation Methods'''
    #--------------------------------------------------------------------------
    def clone_pores(self,pores,apply_label=['clone'],mode='parents'):
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

            - 'parents': (Default) Each clone is connected only to its parent
            - 'siblings': Clones are only connected to each other in the same manner as parents were connected
            - 'isolated': No connections between parents or siblings
        '''
        if (self._geometries != []):
            logger.warning('Network has active Geometries, new pores must be assigned a Geometry')
        if (self._phases != []):
            raise Exception('Network has active Phases, cannot proceed')

        logger.debug('Cloning pores')
        apply_label = list(apply_label)
        #Clone pores
        Np = self.num_pores()
        Nt = self.num_throats()
        parents = sp.array(pores,ndmin=1)
        pcurrent = self['pore.coords']
        pclone = pcurrent[pores,:]
        pnew = sp.concatenate((pcurrent,pclone),axis=0)
        Npnew = sp.shape(pnew)[0]
        clones = sp.arange(Np,Npnew)
        #Add clone labels to network
        for item in apply_label:
            if ('pore.'+item) not in self.keys():
                self['pore.'+item] = False
            if ('throat.'+item) not in self.keys():
                self['throat.'+item] = False
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
        self._update_network()

    def stitch(self,donor,pores_1,pores_2,method='delaunay',len_max=sp.inf,label_suffix=''):
        r'''
        Stitches a second a network to the current network.

        Parameters
        ----------
        donor : OpenPNM Network Object
            The Network to stitch on to the current Network

        pores_1 : array_like
            The pores on the current Network

        pores_2 : array_like
            The pores on the donor Network
            
        label_suffix : string or None
            Some text to append to each label in the donor Network before
            inserting them into the recipient.  The default is to append no 
            text, but a common option would be to append the donor Network's 
            name. To insert none of the donor labels, use None.

        len_max : float
            Set a length limit on length of new throats

        method : string (default = 'delaunay')
            The method to use when making pore to pore connections. Options are:

            - 'delaunay' : Use a Delaunay tessellation
            - 'nearest' : Connects each pore on the receptor network to its nearest pore on the donor network
            
        Notes
        -----
        Before stitching it is necessary to translate the pore coordinates of 
        one of the Networks so that it is positioned correctly relative to the
        other.  
        
        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn2 = OpenPNM.Network.TestNet()
        >>> [pn.Np, pn.Nt]
        [125, 300]
        >>> [pn2.Np, pn2.Nt]
        [125, 300]
        >>> pn2['pore.coords'][:,2] += 5.0  # Translate pn2 up 5 units in the Z-direction
        >>> pn.stitch(donor=pn2,pores_1=pn.pores('top'),pores_2=pn2.pores('bottom'),len_max=1.0)
        >>> [pn.Np, pn.Nt]
        [250, 625]

        '''
        # Ensure Networks have no associated objects yet
        if (len(self._simulation()) > 1) or (len(donor._simulation()) > 1):
            raise Exception('Cannot stitch a Network with active sibling objects')
        # Get the initial number of pores and throats
        N_init = {}
        N_init['pore']  = self.Np
        N_init['throat'] = self.Nt
        if method == 'delaunay':
            P1 = pores_1
            P2 = pores_2 + N_init['pore']  # Increment pores on donor
            P = sp.hstack((P1,P2))
            C1 = self['pore.coords'][pores_1]
            C2 = donor['pore.coords'][pores_2]
            C = sp.vstack((C1,C2))
            T = sp.spatial.Delaunay(C)
            a = T.simplices
            b = []
            for i in range(0,4):  # Turn simplices into pairs of connections
                for j in range(0,4):
                    if i != j:
                        b.append(sp.vstack((a[:,i],a[:,j])).T)
            c = sp.vstack((b))  # Convert list b to numpy array
            #Remove thraot connections to self
            K1 = sp.in1d(P[c][:,0],P1)*sp.in1d(P[c][:,1],P2)
            K2 = sp.in1d(P[c][:,0],P2)*sp.in1d(P[c][:,1],P1)
            K = sp.where(K1 + K2)[0]
            #Remove duplicate throat connections
            i = P[c][K][:,0]
            j = P[c][K][:,1]
            v = sp.ones_like(i)
            N = self.Np + donor.Np
            adjmat = sprs.coo.coo_matrix((v,(i,j)),shape=(N,N))
            #Convert to csr and back to coo to remove duplicates
            adjmat = adjmat.tocsr()
            adjmat = adjmat.tocoo()
            #Remove lower triangular to remove bidirectional
            adjmat = sprs.triu(adjmat,k=1,format="coo")
            #Convert adjmat into 'throat.conns'
            conns = sp.vstack((adjmat.row, adjmat.col)).T
        if method == 'nearest':
            P1 = pores_1
            P2 = pores_2 + N_init['pore'] # Increment pores on donor
            P = sp.hstack((P1,P2))
            C1 = self['pore.coords'][pores_1]
            C2 = donor['pore.coords'][pores_2]
            D = sp.spatial.distance.cdist(C1,C2)
            [P1_ind,P2_ind] = sp.where(D<=len_max)
            conns = sp.vstack((P1[P1_ind],P2[P2_ind])).T

        #Enter donor's pores into the Network
        self.extend(pore_coords=donor['pore.coords'])

        #Enter donor's throats into the Network
        self.extend(throat_conns=donor['throat.conns']+N_init['pore'])

        #Trim throats that are longer then given len_max
        C1 = self['pore.coords'][conns[:,0]]
        C2 = self['pore.coords'][conns[:,1]]
        L = sp.sum((C1 - C2)**2,axis=1)**0.5
#        conns = conns[L<=len_max]

        #Add donor labels to recipient network
        if label_suffix != None:
            if label_suffix != '':
                label_suffix = '_'+label_suffix
            for label in donor.labels():
                element = label.split('.')[0]
                locations = sp.where(self._get_indices(element)>=N_init[element])[0]
                try:
                    self[label+label_suffix]
                except:
                    self[label+label_suffix] = False
                self[label+label_suffix][locations] = donor[label]

        #Add the new stitch throats to the Network
        self.extend(throat_conns=conns,labels='stitched')
        
        # Remove donor from Controller, if present
        # This check allows for the reuse of a donor Network multiple times
        if donor in ctrl.values():
            ctrl.purge_object(donor)

    def check_network_health(self):
        r'''
        This method check the network topological health by checking for:

            (1) Isolated pores
            (2) Islands or isolated clusters of pores
            (3) Duplicate throats
            (4) Bidirectional throats (ie. symmetrical adjacency matrix)

        Returns
        -------
        A dictionary containing the offending pores or throat numbers under
        each named key.

        It also returns a list of which pores and throats should be trimmed
        from the network to restore health.  This list is a suggestion only,
        and is based on keeping the largest cluster and trimming the others.

        Notes
        -----
        - Does not yet check for duplicate pores
        - Does not yet suggest which throats to remove
        - This is just a 'check' method and does not 'fix' the problems it finds
        '''

        health = Tools.HealthDict()
        health['disconnected_clusters'] = []
        health['isolated_pores'] = []
        health['trim_pores'] = []
        health['duplicate_throats'] = []
        health['bidirectional_throats'] = []

        #Check for individual isolated pores
        Ps = self.num_neighbors(self.pores())
        if sp.sum(Ps==0) > 0:
            logger.warning(str(sp.sum(Ps==0))+' pores have no neighbors')
            health['isolated_pores'] = sp.where(Ps==0)[0]

        #Check for separated clusters of pores
        temp = []
        Cs = self.find_clusters(self.tomask(throats=self.throats('all')))
        if sp.shape(sp.unique(Cs))[0] > 1:
            logger.warning('Isolated clusters exist in the network')
            for i in sp.unique(Cs):
                temp.append(sp.where(Cs==i)[0])
            b = sp.array([len(item) for item in temp])
            c = sp.argsort(b)[::-1]
            for i in range(0,len(c)):
                health['disconnected_clusters'].append(temp[c[i]])
                if i > 0:
                    health['trim_pores'].extend(temp[c[i]])

        #Check for duplicate throats
        i = self['throat.conns'][:,0]
        j = self['throat.conns'][:,1]
        v = sp.array(self['throat.all'],dtype=int)
        Np = self.num_pores()
        adjmat = sprs.coo_matrix((v,(i,j)),[Np,Np])
        temp = adjmat.tocsr()  # Convert to CSR to combine duplicates
        temp = adjmat.tocoo()  # And back to COO
        mergedTs = sp.where(temp.data>1)
        Ps12 = sp.vstack((temp.row[mergedTs], temp.col[mergedTs])).T
        dupTs = []
        for i in range(0,sp.shape(Ps12)[0]):
            dupTs.append(self.find_connecting_throat(Ps12[i,0],Ps12[i,1]).tolist)
        health['duplicate_throats'] = dupTs

        #Check for bidirectional throats
        num_full = adjmat.sum()
        temp = sprs.triu(adjmat,k=1)
        num_upper = temp.sum()
        if num_full > num_upper:
            health['bidirectional_throats'] = str(num_full-num_upper)+' detected!'

        #Check for coincident pores
#        temp = misc.dist(self['pore.coords'],self['pore.coords'])
#        temp = sp.triu(temp,k=1)  # Remove lower triangular of matrix
#        temp = sp.where(temp==0)  # Find 0 values in distance matrix
#        dupPs = sp.where(temp[1]>temp[0])[0]  # Find 0 values above diagonal
#        health['duplicate_pores'] = dupPs

        return health

    def check_geometry_health(self):
        r'''
        Perform a check to find pores with overlapping or undefined Geometries
        '''
        geoms = self.geometries()
        Ptemp = sp.zeros((self.Np,))
        Ttemp = sp.zeros((self.Nt,))
        for item in geoms:
            Pind = self['pore.'+item]
            Tind = self['throat.'+item]
            Ptemp[Pind] = Ptemp[Pind] + 1
            Ttemp[Tind] = Ttemp[Tind] + 1
        health = Tools.HealthDict()
        health['overlapping_pores'] = sp.where(Ptemp>1)[0].tolist()
        health['undefined_pores'] = sp.where(Ptemp==0)[0].tolist()
        health['overlapping_throats'] = sp.where(Ttemp>1)[0].tolist()
        health['undefined_throats'] = sp.where(Ttemp==0)[0].tolist()
        return health

    def _update_network(self,mode='clear'):
        r'''
        Regenerates the adjacency and incidence matrices

        Parameters
        ----------
        mode : string
            Controls the extent of the update.  Options are:

            - 'clear' : Removes exsiting adjacency and incidence matrices
            - 'regenerate' : Removes the existing matrices and regenerates new ones.

        Notes
        -----
        The 'regenerate' mode is more time consuming, so repeated calls to
        this function (ie. during network merges, and adding boundaries)
        should use the 'clear' mode.  The other methods that require these
        matrices will generate them as needed, so this pushes the 'generation'
        time to 'on demand'.
        '''
        logger.debug('Resetting adjacency and incidence matrices')
        self._adjacency_matrix['coo'] = {}
        self._adjacency_matrix['csr'] = {}
        self._adjacency_matrix['lil'] = {}
        self._incidence_matrix['coo'] = {}
        self._incidence_matrix['csr'] = {}
        self._incidence_matrix['lil'] = {}

        if mode == 'regenerate':
            self._adjacency_matrix['coo'] = self.create_adjacency_matrix(sprsfmt='coo')
            self._adjacency_matrix['csr'] = self.create_adjacency_matrix(sprsfmt='csr')
            self._adjacency_matrix['lil'] = self.create_adjacency_matrix(sprsfmt='lil')
            self._incidence_matrix['coo'] = self.create_incidence_matrix(sprsfmt='coo')
            self._incidence_matrix['csr'] = self.create_incidence_matrix(sprsfmt='csr')
            self._incidence_matrix['lil'] = self.create_incidence_matrix(sprsfmt='lil')

    #--------------------------------------------------------------------------
    '''Domain Geometry Methods'''
    #--------------------------------------------------------------------------
    def domain_bulk_volume(self):
        r'''
        '''
        raise NotImplementedError()

    def isolated_pores(self):
        r'''
        This method checks to see whether any pores are isolated from the network and
        returns a boolean mask
        '''
        isolated = [False]*(self.num_pores())
        for pore in self.pores():
            if pore not in self["throat.conns"]:
                isolated[pore]=True
        return isolated

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
            logger.warning('The supplied pores are not coplanar. Length will be approximate.')
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
        coords = self['pore.coords'][face]
        rads = self['pore.diameter'][face]/2.
        # calculate the area of the 3 principle faces of the bounding cuboid
        dx = max(coords[:,0]+rads) - min(coords[:,0]-rads)
        dy = max(coords[:,1]+rads) - min(coords[:,1]-rads)
        dz = max(coords[:,2]+rads) - min(coords[:,2]-rads)
        yz = dy*dz # x normal
        xz = dx*dz # y normal
        xy = dx*dy # z normal
        # find the directions parallel to the plane
        directions = sp.where([yz,xz,xy]!=max([yz,xz,xy]))[0]
        try:
            # now, use the whole network to do the area calculation
            coords = self['pore.coords']
            rads = self['pore.diameter']/2.
            d0 = (max(coords[:,directions[0]]+rads) - min(coords[:,directions[0]]-rads))
            d1 = (max(coords[:,directions[1]]+rads) - min(coords[:,directions[1]]-rads))
            A = d0*d1
        except:
            # if that fails, use the max face area of the bounding cuboid
            A = max([yz,xz,xy])
        if not misc.iscoplanar(self['pore.coords'][face]):
            logger.warning('The supplied pores are not coplanar. Area will be approximate')
            pass
        return A

if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)










