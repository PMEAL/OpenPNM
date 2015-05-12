# -*- coding: utf-8 -*-
"""
===============================================================================
GenericNetwork: Abstract class to construct pore networks
===============================================================================

"""
import scipy as sp
import scipy.sparse as sprs
import scipy.spatial as sptl
import OpenPNM.Utilities.misc as misc
from OpenPNM.Utilities import topology
from OpenPNM.Base import Core, Controller, Tools, logging
logger = logging.getLogger(__name__)
ctrl = Controller()
topo = topology()


class GenericNetwork(Core):
    r"""
    GenericNetwork - Base class to construct pore networks

    Parameters
    ----------
    name : string
        Unique name for Network object

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        logger.name = self.name

        # Initialize adjacency and incidence matrix dictionaries
        self._incidence_matrix = {}
        self._adjacency_matrix = {}

    def __setitem__(self, prop, value):
        if prop == 'throat.conns':
            if sp.shape(value)[1] != 2:
                logger.error('Wrong size for throat conns!')
            else:
                mask = value[:, 0] > value[:, 1]
                if mask.any():
                    logger.debug('The first column in (throat.conns) should be \
                                  smaller than the second one.')
                    v1 = sp.copy(value[:, 0][mask])
                    v2 = sp.copy(value[:, 1][mask])
                    value[:, 0][mask] = v2
                    value[:, 1][mask] = v1
        for geom in self._geometries:
            if (prop in geom.keys()) and ('all' not in prop.split('.')):
                logger.error(prop + ' is already defined in at least one associated \
                             Geometry object')
                return
        super().__setitem__(prop, value)

    def __getitem__(self, key):
        if key.split('.')[-1] == self.name:
            element = key.split('.')[0]
            return self[element+'.all']
        if key not in self.keys():
            logger.debug(key + ' not on Network, constructing data from Geometries')
            return self._interleave_data(key, self.geometries())
        else:
            return super().__getitem__(key)

    def _set_net(self, network):
        pass

    def _get_net(self):
        return self
    _net = property(fset=_set_net, fget=_get_net)

    def create_adjacency_matrix(self, data=None, sprsfmt='coo',
                                dropzeros=True, sym=True):
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
        >>> temp = pn.create_adjacency_matrix(data=vals, sprsfmt='csr')

        """
        logger.debug('create_adjacency_matrix: Start of method')
        Np = self.num_pores()
        Nt = self.num_throats()

        # Check if provided data is valid
        if data is None:
            data = sp.ones((self.num_throats(),))
        elif sp.shape(data)[0] != Nt:
            raise Exception('Received dataset of incorrect length')

        # Clear any zero-weighted connections
        if dropzeros:
            ind = data > 0
        else:
            ind = sp.ones_like(data, dtype=bool)

        # Get connectivity info from network
        conn = self['throat.conns'][ind]
        row = conn[:, 0]
        col = conn[:, 1]
        data = data[ind]

        # Append row & col to each other, and data to itself
        if sym:
            row = sp.append(row, conn[:, 1])
            col = sp.append(col, conn[:, 0])
            data = sp.append(data, data)

        # Generate sparse adjacency matrix in 'coo' format
        temp = sprs.coo_matrix((data, (row, col)), (Np, Np))

        # Convert to requested format
        if sprsfmt == 'coo':
            pass  # temp is already in coo format
        if sprsfmt == 'csr':
            temp = temp.tocsr()
        if sprsfmt == 'lil':
            temp = temp.tolil()
        logger.debug('create_adjacency_matrix: End of method')
        return temp

    def create_incidence_matrix(self, data=None, sprsfmt='coo', dropzeros=True):
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

        # Check if provided data is valid
        if data is None:
            data = sp.ones((self.num_throats(),))
        elif sp.shape(data)[0] != Nt:
            raise Exception('Received dataset of incorrect length')

        if dropzeros:
            ind = data > 0
        else:
            ind = sp.ones_like(data, dtype=bool)

        conn = self['throat.conns'][ind]
        row = conn[:, 0]
        row = sp.append(row, conn[:, 1])
        col = self.throats('all')[ind]
        col = sp.append(col, col)
        data = sp.append(data[ind], data[ind])

        temp = sprs.coo.coo_matrix((data, (row, col)), (Np, Nt))

        # Convert to requested format
        if sprsfmt == 'coo':
            pass  # temp is already in coo format
        if sprsfmt == 'csr':
            temp = temp.tocsr()
        if sprsfmt == 'lil':
            temp = temp.tolil()
        logger.debug('create_incidence_matrix: End of method')
        return temp

    def find_connected_pores(self, throats=[], flatten=False):
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
        >>> pn.find_connected_pores(throats=[0,1], flatten=True)
        array([0, 1, 5])
        """
        Ps = self['throat.conns'][throats]
        if flatten:
            Ps = sp.unique(sp.hstack(Ps))
        return Ps

    def find_connecting_throat(self, P1, P2):
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
        >>> pn.find_connecting_throat([0, 1, 2], [2, 2, 2])
        [[], [3], []]

        TODO: This now works on 'vector' inputs, but is not actually vectorized
        in the Numpy sense, so could be slow with large P1,P2 inputs
        """
        Ts1 = self.find_neighbor_throats(P1, flatten=False)
        Ts2 = self.find_neighbor_throats(P2, flatten=False)
        Ts = []
        if sp.shape(P1) == ():
            P1 = [P1]
            P2 = [P2]
        for row in range(0, len(P1)):
            if P1[row] == P2[row]:
                throat = []
            else:
                throat = sp.intersect1d(Ts1[row], Ts2[row]).tolist()
            Ts.insert(0, throat)
        Ts.reverse()
        return Ts

    def find_neighbor_pores(self, pores, mode='union', flatten=True, excl_self=True):
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
        >>> pn.find_neighbor_pores(pores=[0, 2])
        array([ 1,  3,  5,  7, 25, 27])
        >>> pn.find_neighbor_pores(pores=[0, 1])
        array([ 2,  5,  6, 25, 26])
        >>> pn.find_neighbor_pores(pores=[0, 1], mode='union', excl_self=False)
        array([ 0,  1,  2,  5,  6, 25, 26])
        >>> pn.find_neighbor_pores(pores=[0, 2],flatten=False)
        array([array([ 1,  5, 25]), array([ 1,  3,  7, 27])], dtype=object)
        >>> pn.find_neighbor_pores(pores=[0, 2],mode='intersection')
        array([1])
        >>> pn.find_neighbor_pores(pores=[0, 2],mode='not_intersection')
        array([ 3,  5,  7, 25, 27])
        """
        pores = sp.array(pores, ndmin=1)
        try:
            neighborPs = self._adjacency_matrix['lil'].rows[[pores]]
        except:
            temp = self.create_adjacency_matrix(sprsfmt='lil')
            self._adjacency_matrix['lil'] = temp
            neighborPs = self._adjacency_matrix['lil'].rows[[pores]]
        if [sp.asarray(x) for x in neighborPs if x] == []:
            return sp.array([], ndmin=1)
        if flatten:
            # All the empty lists must be removed to maintain data type after
            # hstack (numpy bug?)
            neighborPs = [sp.asarray(x) for x in neighborPs if x]
            neighborPs = sp.hstack(neighborPs)
            neighborPs = sp.concatenate((neighborPs, pores))
            # Remove references to input pores and duplicates
            if mode == 'not_intersection':
                neighborPs = sp.array(sp.unique(sp.where(
                    sp.bincount(neighborPs) == 1)[0]), dtype=int)
            elif mode == 'union':
                neighborPs = sp.array(sp.unique(neighborPs), int)
            elif mode == 'intersection':
                neighborPs = sp.array(sp.unique(sp.where(
                    sp.bincount(neighborPs) > 1)[0]), dtype=int)
            if excl_self:
                neighborPs = neighborPs[~sp.in1d(neighborPs, pores)]
        else:
            for i in range(0, sp.size(pores)):
                neighborPs[i] = sp.array(neighborPs[i], dtype=int)
        return sp.array(neighborPs, ndmin=1)

    def find_neighbor_throats(self, pores, mode='union', flatten=True):
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
        >>> pn.find_neighbor_throats(pores=[0, 1])
        array([0, 1, 2, 3, 4, 5])
        >>> pn.find_neighbor_throats(pores=[0, 1],flatten=False)
        array([array([0, 1, 2]), array([0, 3, 4, 5])], dtype=object)
        """
        # Test for existence of incidence matrix
        try:
            neighborTs = self._incidence_matrix['lil'].rows[[pores]]
        except:
            temp = self.create_incidence_matrix(sprsfmt='lil')
            self._incidence_matrix['lil'] = temp
            neighborTs = self._incidence_matrix['lil'].rows[[pores]]
        if [sp.asarray(x) for x in neighborTs if x] == []:
            return sp.array([], ndmin=1)
        if flatten:
            # All the empty lists must be removed to maintain data type after
            # hstack (numpy bug?)
            neighborTs = [sp.asarray(x) for x in neighborTs if x]
            neighborTs = sp.hstack(neighborTs)
            # Remove references to input pores and duplicates
            if mode == 'not_intersection':
                neighborTs = sp.unique(sp.where(sp.bincount(neighborTs) == 1)[0])
            elif mode == 'union':
                neighborTs = sp.unique(neighborTs)
            elif mode == 'intersection':
                neighborTs = sp.unique(sp.where(sp.bincount(neighborTs) > 1)[0])
        else:
            for i in range(0, sp.size(pores)):
                neighborTs[i] = sp.array(neighborTs[i])
        return sp.array(neighborTs, ndmin=1)

    def num_neighbors(self, pores, flatten=False):
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
        >>> pn.num_neighbors(pores=[0, 1], flatten=False)
        array([3, 4])
        >>> pn.num_neighbors(pores=[0, 1], flatten=True)
        5
        >>> pn.num_neighbors(pores=[0, 2], flatten=True)
        6
        """

        # Count number of neighbors
        if flatten:
            neighborPs = self.find_neighbor_pores(pores, flatten=True,
                                                  mode='union', excl_self=True)
            num = sp.shape(neighborPs)[0]
        else:
            neighborPs = self.find_neighbor_pores(pores, flatten=False)
            num = sp.zeros(sp.shape(neighborPs), dtype=int)
            for i in range(0, sp.shape(num)[0]):
                num[i] = sp.size(neighborPs[i])
        return num

    def find_interface_throats(self, labels=[]):
        r"""
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
        >>> pn['pore.domain1'][[0, 1, 2]] = True
        >>> pn['pore.domain2'][[5, 6, 7]] = True
        >>> pn.find_interface_throats(labels=['domain1', 'domain2'])
        array([1, 4, 7])
        """
        Tind = sp.array([], ndmin=1)
        if sp.shape(labels)[0] != 2:
            logger.error('Exactly two labels must be given')
            pass
        else:
            P1 = self.pores(labels=labels[0])
            P2 = self.pores(labels=labels[1])
            # Check if labels overlap
            if sp.sum(sp.in1d(P1, P2)) > 0:
                logger.error('Some labels overlap, iterface cannot be found')
                pass
            else:
                T1 = self.find_neighbor_throats(P1)
                T2 = self.find_neighbor_throats(P2)
                Tmask = sp.in1d(T1, T2)
                Tind = T1[Tmask]
        return Tind

    def find_clusters(self, mask=[]):
        r"""
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

        """
        if sp.shape(mask)[0] == self.num_throats():
            # Convert to boolean mask if not already
            temp = sp.zeros((self.num_throats(),), dtype=bool)
            temp[mask] = True
        elif sp.shape(mask)[0] == self.num_pores():
            conns = self.find_connected_pores(throats=self.throats())
            conns[:, 0] = mask[conns[:, 0]]
            conns[:, 1] = mask[conns[:, 1]]
            temp = sp.array(conns[:, 0]*conns[:, 1], dtype=bool)
        else:
            raise Exception('Mask received was neither Nt nor Np long')
        temp = self.create_adjacency_matrix(data=temp, sprsfmt='csr', dropzeros=True)
        clusters = sprs.csgraph.connected_components(csgraph=temp,
                                                     directed=False)[1]
        return clusters

    def find_nearest_pores(self, pores, distance=0):
        r"""
        Find all pores with a given distance of the input pore(s) regardless of
        whether or not they are toplogically connected.

        Parameters
        ----------
        pores : array_like
            The list of pores for whom nearby neighbors are to be found
        distance : scalar
            The within which the nearby should be found

        Returns
        -------
            A list of pores which are within the given spatial distance.  If a
            list of N pores is supplied, then a an N-long list of such lists is
            returned.  The returned lists each contain the pore for which the
            neighbors were sought.

        Examples
        --------
        >>> import OpenPNM as op
        >>> pn = op.Network.TestNet()
        >>> pn.find_nearest_pores(pores=[0, 1],distance=1)
        array([[0, 1, 25, 5], [0, 1, 26, 6, 2]], dtype=object)
        >>> pn.find_nearest_pores(pores=[0, 1],distance=.5)
        array([[0], [1]], dtype=object)
        """
        kd = sptl.cKDTree(self['pore.coords'])
        if distance == 0:
            pass
        elif distance > 0:
            Pn = kd.query_ball_point(self['pore.coords'][pores], r=distance)
        return Pn

    def extend(self, pore_coords=[], throat_conns=[], labels=[]):
        topo.extend(network=self, pore_coords=pore_coords,
                    throat_conns=throat_conns, labels=labels)
    extend.__doc__ = topo.extend.__doc__

    def trim(self, pores=[], throats=[]):
        topo.trim(network=self, pores=pores, throats=throats)
    trim.__doc__ = topo.trim.__doc__

    def clone_pores(self, pores, apply_label=['clone'], mode='parents'):
        topo.clone_pores(network=self, pores=pores,
                         apply_label=apply_label, mode=mode)
    clone_pores.__doc__ = topo.clone_pores.__doc__

    def stitch(self, donor, P_donor, P_network, len_max=sp.inf, label_suffix=''):
        topo.stitch(network=self, donor=donor, P_donor=P_donor, P_network=P_network,
                    len_max=len_max, label_suffix=label_suffix)
    stitch.__doc__ = topo.stitch.__doc__

    def connect_pores(self, pores1, pores2, labels=[]):
        topo.connect_pores(network=self, pores1=pores1, pores2=pores2, labels=labels)
    connect_pores.__doc__ = topo.connect_pores.__doc__

    def check_network_health(self):
        r"""
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
        """

        health = Tools.HealthDict()
        health['disconnected_clusters'] = []
        health['isolated_pores'] = []
        health['trim_pores'] = []
        health['duplicate_throats'] = []
        health['bidirectional_throats'] = []

        # Check for individual isolated pores
        Ps = self.num_neighbors(self.pores())
        if sp.sum(Ps == 0) > 0:
            logger.warning(str(sp.sum(Ps == 0)) + ' pores have no neighbors')
            health['isolated_pores'] = sp.where(Ps == 0)[0]

        # Check for separated clusters of pores
        temp = []
        Cs = self.find_clusters(self.tomask(throats=self.throats('all')))
        if sp.shape(sp.unique(Cs))[0] > 1:
            logger.warning('Isolated clusters exist in the network')
            for i in sp.unique(Cs):
                temp.append(sp.where(Cs == i)[0])
            b = sp.array([len(item) for item in temp])
            c = sp.argsort(b)[::-1]
            for i in range(0, len(c)):
                health['disconnected_clusters'].append(temp[c[i]])
                if i > 0:
                    health['trim_pores'].extend(temp[c[i]])

        # Check for duplicate throats
        i = self['throat.conns'][:, 0]
        j = self['throat.conns'][:, 1]
        v = sp.array(self['throat.all'], dtype=int)
        Np = self.num_pores()
        adjmat = sprs.coo_matrix((v, (i, j)), [Np, Np])
        temp = adjmat.tocsr()  # Convert to CSR to combine duplicates
        temp = adjmat.tocoo()  # And back to COO
        mergedTs = sp.where(temp.data > 1)
        Ps12 = sp.vstack((temp.row[mergedTs], temp.col[mergedTs])).T
        dupTs = []
        for i in range(0, sp.shape(Ps12)[0]):
            dupTs.append(self.find_connecting_throat(Ps12[i, 0], Ps12[i, 1]).tolist)
        health['duplicate_throats'] = dupTs

        # Check for bidirectional throats
        num_full = adjmat.sum()
        temp = sprs.triu(adjmat, k=1)
        num_upper = temp.sum()
        if num_full > num_upper:
            health['bidirectional_throats'] = str(num_full-num_upper) + ' detected!'

        return health

    def check_geometry_health(self):
        r"""
        Perform a check to find pores with overlapping or undefined Geometries
        """
        geoms = self.geometries()
        Ptemp = sp.zeros((self.Np,))
        Ttemp = sp.zeros((self.Nt,))
        for item in geoms:
            Pind = self['pore.'+item]
            Tind = self['throat.'+item]
            Ptemp[Pind] = Ptemp[Pind] + 1
            Ttemp[Tind] = Ttemp[Tind] + 1
        health = Tools.HealthDict()
        health['overlapping_pores'] = sp.where(Ptemp > 1)[0].tolist()
        health['undefined_pores'] = sp.where(Ptemp == 0)[0].tolist()
        health['overlapping_throats'] = sp.where(Ttemp > 1)[0].tolist()
        health['undefined_throats'] = sp.where(Ttemp == 0)[0].tolist()
        return health

    def _update_network(self, mode='clear'):
        r"""
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
        """
        logger.debug('Resetting adjacency and incidence matrices')
        self._adjacency_matrix['coo'] = {}
        self._adjacency_matrix['csr'] = {}
        self._adjacency_matrix['lil'] = {}
        self._incidence_matrix['coo'] = {}
        self._incidence_matrix['csr'] = {}
        self._incidence_matrix['lil'] = {}

        if mode == 'regenerate':
            self._adjacency_matrix['coo'] = \
                self.create_adjacency_matrix(sprsfmt='coo')
            self._adjacency_matrix['csr'] = \
                self.create_adjacency_matrix(sprsfmt='csr')
            self._adjacency_matrix['lil'] = \
                self.create_adjacency_matrix(sprsfmt='lil')
            self._incidence_matrix['coo'] = \
                self.create_incidence_matrix(sprsfmt='coo')
            self._incidence_matrix['csr'] = \
                self.create_incidence_matrix(sprsfmt='csr')
            self._incidence_matrix['lil'] = \
                self.create_incidence_matrix(sprsfmt='lil')

    def domain_bulk_volume(self):
        raise NotImplementedError()

    def domain_pore_volume(self):
        raise NotImplementedError()

    def domain_length(self, face_1, face_2):
        r"""
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
        """
        # Ensure given points are coplanar before proceeding
        if misc.iscoplanar(self['pore.coords'][face_1]) and \
                misc.iscoplanar(self['pore.coords'][face_2]):
            # Find distance between given faces
            x = self['pore.coords'][face_1]
            y = self['pore.coords'][face_2]
            Ds = misc.dist(x, y)
            L = sp.median(sp.amin(Ds, axis=0))
        else:
            logger.warning('The supplied pores are not coplanar. Length will be \
                            approximate.')
            f1 = self['pore.coords'][face_1]
            f2 = self['pore.coords'][face_2]
            distavg = [0, 0, 0]
            distavg[0] = sp.absolute(sp.average(f1[:, 0]) - sp.average(f2[:, 0]))
            distavg[1] = sp.absolute(sp.average(f1[:, 1]) - sp.average(f2[:, 1]))
            distavg[2] = sp.absolute(sp.average(f1[:, 2]) - sp.average(f2[:, 2]))
            L = max(distavg)
        return L

    def domain_area(self, face):
        r"""
        Calculate the area of a given network face

        Parameters
        ----------
        face : array_like
            List of pores of pore defining the face of interest

        Returns
        -------
        The area of the specified face
        """
        coords = self['pore.coords'][face]
        rads = self['pore.diameter'][face] / 2.
        # Calculate the area of the 3 principle faces of the bounding cuboid
        dx = max(coords[:, 0]+rads) - min(coords[:, 0] - rads)
        dy = max(coords[:, 1]+rads) - min(coords[:, 1] - rads)
        dz = max(coords[:, 2]+rads) - min(coords[:, 2] - rads)
        yz = dy*dz  # x normal
        xz = dx*dz  # y normal
        xy = dx*dy  # z normal
        # Find the directions parallel to the plane
        directions = sp.where([yz, xz, xy] != max([yz, xz, xy]))[0]
        try:
            # Use the whole network to do the area calculation
            coords = self['pore.coords']
            rads = self['pore.diameter']/2.
            d0 = max(coords[:, directions[0]] + rads) - \
                min(coords[:, directions[0]] - rads)
            d1 = max(coords[:, directions[1]] + rads) - \
                min(coords[:, directions[1]] - rads)
            A = d0*d1
        except:
            # If that fails, use the max face area of the bounding cuboid
            A = max([yz, xz, xy])
        if not misc.iscoplanar(self['pore.coords'][face]):
            logger.warning('The supplied pores are not coplanar. Area will be \
                            approximate')
            pass
        return A
