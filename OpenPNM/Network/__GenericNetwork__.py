# -*- coding: utf-8 -*-
"""
===============================================================================
GenericNetwork: Abstract class to construct pore networks
===============================================================================

"""
import itertools
import scipy as sp
import scipy.sparse as sprs
import scipy.spatial as sptl
import OpenPNM.Utilities.misc as misc
from OpenPNM.Utilities import topology
from OpenPNM.Base import Core, Workspace, Tools, logging
logger = logging.getLogger(__name__)
mgr = Workspace()
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
        self.network.update({self.name: self})

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

    def props(self, element=None, mode='all', deep=False):
        # TODO: The mode 'deep' is deprecated in favor of the deep argument
        # and should be removed in a future version
        modes = ['all', 'deep', 'models', 'constants']
        mode = self._parse_mode(mode=mode, allowed=modes, single=False)
        prop_list = Tools.PrintableList()
        if ('deep' in mode) or (deep is True):
            if mode.count('deep') > 0:
                mode.remove('deep')
            for geom in self._geometries:
                prop_list.extend(geom.props(element=element, mode=mode))
            # Get unique values
            prop_list = Tools.PrintableList(set(prop_list))
        prop_list.extend(super().props(element=element, mode=mode))
        return prop_list

    props.__doc__ = Core.props.__doc__

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

            **'coo'** : (default) This is the native format of OpenPNM data

            **'lil'** : Enables row-wise slice of data

            **'csr'** : Favored by most linear algebra routines

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

        # Check if provided data is valid
        if data is None:
            data = sp.ones((self.Nt,))
        elif sp.shape(data)[0] != self.Nt:
            raise Exception('Received dataset of incorrect length')

        # Clear any zero-weighted connections
        if dropzeros:
            ind = data > 0
            conn = self['throat.conns'][ind]
            data = data[ind]
        else:
            conn = self['throat.conns']

        # Get connectivity info from network
        row = conn[:, 0]
        col = conn[:, 1]

        # Append row & col to each other, and data to itself
        if sym:
            row = sp.append(row, conn[:, 1])
            col = sp.append(col, conn[:, 0])
            data = sp.append(data, data)

        # Generate sparse adjacency matrix in 'coo' format
        temp = sprs.coo_matrix((data, (row, col)), (self.Np, self.Np))

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

            **'coo'** : (default) This is the native format of OpenPNMs data

            **'lil'** : Enables row-wise slice of data

            **'csr'** : Favored by most linear algebra routines

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

        # Check if provided data is valid
        if data is None:
            data = sp.ones((self.Nt,))
        elif sp.shape(data)[0] != self.Nt:
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

        temp = sprs.coo.coo_matrix((data, (row, col)), (self.Np, self.Nt))

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

        Notes
        -----
        This method basically just looks into the pn['throat.conns'] array and
        retrieves the pores for each input throat.  The flatten option merely
        stacks the two columns and eliminate non-unique values.
        """
        Ts = self._parse_locations(throats)
        Ps = self['throat.conns'][Ts]
        if flatten:
            if sp.shape(Ps) == (0, 2):
                Ps = sp.array([], ndmin=1, dtype=int)
            else:
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
        P1 = self._parse_locations(P1)
        P2 = self._parse_locations(P2)
        Ts1 = self.find_neighbor_throats(P1, flatten=False)
        Ts2 = self.find_neighbor_throats(P2, flatten=False)
        Ts = []

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
        excl_self : bool
            If this is True (default) then the input pores are not included in
            the returned list.  This option only applies when input pores are
            in fact neighbors to each other, otherwise they are not part of the
            returned list anyway.
        mode : string, optional
            Specifies which neighbors should be returned.  The options are:

            **'union'** : All neighbors of the input pores

            **'intersection'** : Only neighbors shared by all input pores

            **'not_intersection'** : Only neighbors not shared by any input
            pores

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
        >>> pn.find_neighbor_pores(pores=[0, 2], flatten=False)
        array([array([ 1,  5, 25]), array([ 1,  3,  7, 27])], dtype=object)
        >>> pn.find_neighbor_pores(pores=[0, 2], mode='intersection')
        array([1])
        >>> pn.find_neighbor_pores(pores=[0, 2], mode='not_intersection')
        array([ 3,  5,  7, 25, 27])
        """
        neighbors = self._find_neighbors(pores=pores, element='pore',
                                         mode=mode, flatten=flatten,
                                         excl_self=excl_self)
        return neighbors

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

            **'union'** : All neighbors of the input pores

            **'intersection'** : Only neighbors shared by all input pores

            **'not_intersection'** : Only neighbors not shared by any input
            pores

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
        neighbors = self._find_neighbors(pores=pores, mode=mode,
                                         element='throat', flatten=flatten,
                                         excl_self=False)
        return neighbors

    def _find_neighbors(self, pores, element, mode, flatten, excl_self):
        r"""
        Private method for finding the neighboring pores or throats connected
        directly to given set of pores.

        Parameters
        ----------
        pores : array_like
            The list of pores whose neighbors are sought
        element : string, either 'pore' or 'throat'
            Whether to find neighboring pores or throats
        mode : string
            Controls how the neighbors are filtered.  Options are:

            **'union'** : All neighbors of the input pores

            **'intersection'** : Only neighbors shared by all input pores

            **'not_intersection'** : Only neighbors not shared by any input
            pores

        flatten : boolean
            If flatten is True (default) a 1D array of unique neighbors is
            returned. If flatten is False the returned array contains arrays
            of neighboring throat ID numbers for each input pore, in the order
            they were sent.
        excl_self : bool
            When True the input pores are not included in the returned list of
            neighboring pores.  This option only applies when input pores are
            in fact neighbors to each other, otherwise they are not part of the
            returned list anyway.  This is ignored with the element is
            'throats'.

        See Also
        --------
        find_neighbor_pores
        find_neighbor_throats
        num_neighors

        """
        element = self._parse_element(element=element, single=True)
        pores = self._parse_locations(pores)
        if sp.size(pores) == 0:
            return sp.array([], ndmin=1, dtype=int)

        # Test for existence of incidence or adjacency matrix
        if element == 'pore':
            try:
                neighbors = self._adjacency_matrix['lil'].rows[[pores]]
            except:
                temp = self.create_adjacency_matrix(sprsfmt='lil')
                self._adjacency_matrix['lil'] = temp
                neighbors = self._adjacency_matrix['lil'].rows[[pores]]
        elif element == 'throat':
            try:
                neighbors = self._incidence_matrix['lil'].rows[[pores]]
            except:
                temp = self.create_incidence_matrix(sprsfmt='lil')
                self._incidence_matrix['lil'] = temp
                neighbors = self._incidence_matrix['lil'].rows[[pores]]

        if flatten:
            # Convert rows of lil into single flat list
            neighbors = itertools.chain.from_iterable(neighbors)
            if element == 'pore':  # Add input pores to list
                neighbors = itertools.chain.from_iterable([neighbors, pores])
            # Convert list to numpy array
            neighbors = sp.fromiter(neighbors, dtype=int)
            if mode == 'not_intersection':
                neighbors = sp.unique(sp.where(sp.bincount(neighbors) == 1)[0])
            elif mode == 'union':
                neighbors = sp.unique(neighbors)
            elif mode == 'intersection':
                neighbors = sp.unique(sp.where(sp.bincount(neighbors) > 1)[0])
            if excl_self and element == 'pore':  # Remove input pores from list
                neighbors = neighbors[~sp.in1d(neighbors, pores)]
            return sp.array(neighbors, ndmin=1, dtype=int)
        else:
            # Convert lists in array to numpy arrays
            neighbors = [sp.array(neighbors[i]) for i in range(0, len(pores))]
            return sp.array(neighbors, ndmin=1)

    def num_neighbors(self, pores, element='pore', flatten=False,
                      mode='union'):
        r"""
        Returns an array containing the number of neigbhoring pores or throats
        for each given input pore

        Parameters
        ----------
        pores : array_like
            Pores whose neighbors are to be counted
        flatten : boolean (optional)
            If ``False`` (default) the number of pores neighboring each input
            pore as an array the same length as ``pores``.  If ``True`` the sum
            total number of is counted.
        element : string
            Indicates whether to count number of neighboring pores or throats.
            For some complex networks, such as extracted networks, several
            throats may exist between two pores, so this query will return
            different results depending on whether 'pores' (default) or
            'throats' is specified.
        mode : string (This is ignored if ``flatten`` is False)
            The logic to apply to the returned count of pores.

            **'union'** : (Default) The sum of all neighbors connected all
            input pores.

            **'intersection'** : The number of neighboring pores that are
            shared by all input pores.

            **'not_intersection'** : The number of neighboring pores that are
            NOT shared by any input pores.

        Returns
        -------
        If ``flatten`` is False, a 1D array with number of neighbors in each
        element, otherwise a scalar value of the number of neighbors.

        Notes
        -----
        This method literally just counts the number of elements in the array
        returned by ``find_neighbor_pores`` or ``find_neighbor_throats`` and
        uses the same logic.  Explore those methods if uncertain about the
        meaning of the ``mode`` argument here.

        See Also
        --------
        find_neighbor_pores
        find_neighbor_throats

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
        >>> pn.num_neighbors(pores=[0, 1], element='throat', mode='union',
        ...                  flatten=True)
        6
        """
        pores = self._parse_locations(pores)
        # Count number of neighbors
        num = self._find_neighbors(pores, element=element, flatten=flatten,
                                   mode=mode, excl_self=True)
        num = sp.array([sp.size(i) for i in num], dtype=int)
        if flatten:
            num = sp.sum(num)
            num = int(num)
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

        TODO: It might be a good idea to allow overlapping regions
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

        See Also
        --------
        find_clusters2

        """
        if sp.size(mask) == self.num_throats():
            # Convert to boolean mask if not already
            temp = sp.zeros((self.num_throats(),), dtype=bool)
            temp[mask] = True
        elif sp.size(mask) == self.num_pores():
            conns = self.find_connected_pores(throats=self.throats())
            conns[:, 0] = mask[conns[:, 0]]
            conns[:, 1] = mask[conns[:, 1]]
            temp = sp.array(conns[:, 0]*conns[:, 1], dtype=bool)
        else:
            raise Exception('Mask received was neither Nt nor Np long')
        temp = self.create_adjacency_matrix(data=temp,
                                            sprsfmt='csr',
                                            dropzeros=True)
        clusters = sprs.csgraph.connected_components(csgraph=temp,
                                                     directed=False)[1]
        return clusters

    def find_clusters2(self, mask=[], t_labels=False):
        r"""
        Identify connected clusters of pores in the network.  This method can
        also return a list of throat cluster numbers, which correspond to the
        cluster numbers of the pores to which the throat is connected.  Either
        site and bond percolation can be consider, see description of input
        arguments for details.

        Parameters
        ----------
        mask : array_like, boolean
            A list of active bonds or sites (throats or pores).  If the mask is
            Np long, then the method will perform a site percolation, and if
            the mask is Nt long bond percolation will be performed.

        t_labels : boolean (default id False)
            Indicates if throat cluster numbers should also be returned. If
            true then a tuple containing both p_clusters and t_clusters is
            returned.

        Returns
        -------
        A Np long list of pore clusters numbers, unless t_labels is True in
        which case a tuple containing both pore and throat cluster labels is
        returned.  The label numbers correspond such that pores and throats
        with the same label are part of the same cluster.

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.Cubic(shape=[25, 25, 1])
        >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,
        ...                                         pores=pn.Ps,
        ...                                         throats=pn.Ts)
        >>> geom['pore.seed'] = sp.rand(pn.Np)
        >>> geom['throat.seed'] = sp.rand(pn.Nt)

        Bond percolation is achieved by sending a list of invaded throats:

        >>> (p_bond,t_bond) = pn.find_clusters2(mask=geom['throat.seed'] < 0.3,
        ...                                     t_labels=True)

        Site percolation is achieved by sending a list of invaded pores:

        >>> (p_site,t_site) = pn.find_clusters2(mask=geom['pore.seed'] < 0.3,
        ...                                     t_labels=True)

        To visualize the invasion pattern, use matplotlib's matshow method
        along with the Cubic Network's asarray method which converts list based
        data to square arrays:

        .. code-block:: python

            import matplotlib.pyplot as plt
            im_bond = pn.asarray(p_bond)[:, :, 0]
            im_site = pn.asarray(p_site)[:, :, 0]
            plt.subplot(1, 2, 1)
            plt.imshow(im_site, interpolation='none')
            plt.subplot(1, 2, 2)
            plt.imshow(im_bond, interpolation='none')

        """
        # Parse the input arguments
        mask = sp.array(mask, ndmin=1)
        if mask.dtype != bool:
            raise Exception('Mask must be a boolean array of Np or Nt length')

        # If pore mask was given perform site percolation
        if sp.size(mask) == self.Np:
            (p_clusters, t_clusters) = self._site_percolation(mask)
        # If pore mask was given perform bond percolation
        elif sp.size(mask) == self.Nt:
            (p_clusters, t_clusters) = self._bond_percolation(mask)
        else:
            raise Exception('Mask received was neither Nt nor Np long')

        if t_labels:
            return (p_clusters, t_clusters)
        else:
            return p_clusters

    def _site_percolation(self, pmask):
        r"""
        This private method is called by 'find_clusters2'
        """
        # Find throats that produce site percolation
        conns = sp.copy(self['throat.conns'])
        conns[:, 0] = pmask[conns[:, 0]]
        conns[:, 1] = pmask[conns[:, 1]]
        # Only if both pores are True is the throat set to True
        tmask = sp.all(conns, axis=1)

        # Perform the clustering using scipy.csgraph
        csr = self.create_adjacency_matrix(data=tmask,
                                           sprsfmt='csr',
                                           dropzeros=True)
        clusters = sprs.csgraph.connected_components(csgraph=csr,
                                                     directed=False)[1]

        # Adjust cluster numbers such that non-invaded pores are labelled -1
        # Note: The following line also takes care of assigning cluster numbers
        # to single isolated invaded pores
        p_clusters = (clusters + 1)*(pmask) - 1
        # Label invaded throats with their neighboring pore's label
        t_clusters = clusters[self['throat.conns']]
        ind = (t_clusters[:, 0] == t_clusters[:, 1])
        t_clusters = t_clusters[:, 0]
        # Label non-invaded throats with -1
        t_clusters[~ind] = -1

        return (p_clusters, t_clusters)

    def _bond_percolation(self, tmask):
        r"""
        This private method is called by 'find_clusters2'
        """
        # Perform the clustering using scipy.csgraph
        csr = self.create_adjacency_matrix(data=tmask,
                                           sprsfmt='csr',
                                           dropzeros=True)
        clusters = sprs.csgraph.connected_components(csgraph=csr,
                                                     directed=False)[1]

        # Convert clusters to a more usable output:
        # Find pores attached to each invaded throats
        Ps = self.find_connected_pores(throats=tmask, flatten=True)
        # Adjust cluster numbers such that non-invaded pores are labelled -0
        p_clusters = (clusters + 1)*(self.tomask(pores=Ps).astype(int)) - 1
        # Label invaded throats with their neighboring pore's label
        t_clusters = clusters[self['throat.conns']][:, 0]
        # Label non-invaded throats with -1
        t_clusters[~tmask] = -1

        return (p_clusters, t_clusters)

    def find_nearby_pores(self, pores, distance, flatten=False, excl_self=True):
        r"""
        Find all pores within a given radial distance of the input pore(s)
        regardless of whether or not they are toplogically connected.

        Parameters
        ----------
        pores : array_like
            The list of pores for whom nearby neighbors are to be found
        distance : scalar
            The maximum distance within which the search should be performed
        excl_self : bool
            Controls whether the input pores should be included in the returned
            list.  The default is True which means they are not included.
        flatten : bool
            If true returns a single list of all pores that match the criteria,
            otherwise returns an array containing a sub-array for each input
            pore, where each sub-array contains the pores that are nearby to
            each given input pore.  The default is False.

        Returns
        -------
            A list of pores which are within the given spatial distance.  If a
            list of N pores is supplied, then a an N-long list of such lists is
            returned.  The returned lists each contain the pore for which the
            neighbors were sought.

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn.find_nearby_pores(pores=[0, 1], distance=1)
        array([array([ 1,  5, 25]), array([ 0,  2,  6, 26])], dtype=object)
        >>> pn.find_nearby_pores(pores=[0, 1], distance=0.5)
        array([], shape=(2, 0), dtype=int64)
        """
        pores = self._parse_locations(pores)
        # Handle an empty array if given
        if sp.size(pores) == 0:
            return sp.array([], dtype=sp.int64)
        if distance <= 0:
            logger.error('Provided distances should be greater than 0')
            if flatten:
                Pn = sp.array([])
            else:
                Pn = sp.array([sp.array([]) for i in range(0, len(pores))])
            return Pn.astype(sp.int64)
        # Create kdTree objects
        kd = sptl.cKDTree(self['pore.coords'],
                          balanced_tree=False,
                          compact_nodes=False)
        kd_pores = sptl.cKDTree(self['pore.coords'][pores],
                                balanced_tree=False,
                                compact_nodes=False)
        # Perform search
        Pn = kd_pores.query_ball_tree(kd, r=distance)
        # Sort the indices in each list
        [Pn[i].sort() for i in range(0, sp.size(pores))]
        if flatten:  # Convert list of lists to a flat nd-array
            temp = []
            [temp.extend(Ps) for Ps in Pn]
            Pn = sp.unique(temp)
            if excl_self:  # Remove inputs if necessary
                Pn = Pn[~sp.in1d(Pn, pores)]
        else:  # Convert list of lists to an nd-array of nd-arrays
            if excl_self:  # Remove inputs if necessary
                [Pn[i].remove(pores[i]) for i in range(0, sp.size(pores))]
            temp = []
            [temp.append(sp.array(Pn[i])) for i in range(0, sp.size(pores))]
            Pn = sp.array(temp)
        if Pn.dtype == float:
            Pn = Pn.astype(sp.int64)
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

    def stitch(self, donor, P_donor, P_network, method, len_max=sp.inf,
               label_suffix=''):
        topo.stitch(network=self, donor=donor, P_donor=P_donor,
                    P_network=P_network, method=method, len_max=len_max,
                    label_suffix=label_suffix)
    stitch.__doc__ = topo.stitch.__doc__

    def connect_pores(self, pores1, pores2, labels=[]):
        topo.connect_pores(network=self,
                           pores1=pores1,
                           pores2=pores2,
                           labels=labels)
    connect_pores.__doc__ = topo.connect_pores.__doc__

    def check_network_health(self):
        r"""
        This method check the network topological health by checking for:

            (1) Isolated pores
            (2) Islands or isolated clusters of pores
            (3) Duplicate throats
            (4) Bidirectional throats (ie. symmetrical adjacency matrix)
            (5) Headless throats

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
        health['headless_throats'] = []
        health['looped_throats'] = []

        # Check for headless throats
        hits = sp.where(self['throat.conns'] > self.Np - 1)[0]
        if sp.size(hits) > 0:
            health['headless_throats'] = sp.unique(hits)
            logger.warning('Health check cannot complete due to connectivity '
                           'errors. Please correct existing errors & recheck.')
            return health

        # Check for throats that loop back onto the same pore
        P12 = self['throat.conns']
        hits = sp.where(P12[:, 0] == P12[:, 1])[0]
        if sp.size(hits) > 0:
            health['looped_throats'] = hits

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
        am = self.create_adjacency_matrix(sprsfmt='lil')
        mergeTs = []
        for i in range(0, self.Np):
            for j in sp.where(sp.array(am.data[i]) > 1)[0]:
                k = am.rows[i][j]
                mergeTs.extend(self.find_connecting_throat(i, k))
        # Remove duplicates
        mergeTs = [list(i) for i in set(tuple(i) for i in mergeTs)]
        health['duplicate_throats'] = mergeTs

        # Check for bidirectional throats
        adjmat = self.create_adjacency_matrix(sprsfmt='coo')
        num_full = adjmat.sum()
        temp = sprs.triu(adjmat, k=1)
        num_upper = temp.sum()
        if num_full > num_upper:
            biTs = sp.where(self['throat.conns'][:, 0] >
                            self['throat.conns'][:, 1])[0]
            health['bidirectional_throats'] = biTs.tolist()

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
            self._adjacency_matrix['csr'] = self._adjacency_matrix['coo'].tocsr()
            self._adjacency_matrix['lil'] = self._adjacency_matrix['coo'].tolil()
            self._incidence_matrix['coo'] = \
                self.create_incidence_matrix(sprsfmt='coo')
            self._incidence_matrix['csr'] = self._incidence_matrix['coo'].tocsr()
            self._incidence_matrix['lil'] = self._incidence_matrix['coo'].tolil()

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
            distavg[0] = sp.absolute(sp.average(f1[:, 0])-sp.average(f2[:, 0]))
            distavg[1] = sp.absolute(sp.average(f1[:, 1])-sp.average(f2[:, 1]))
            distavg[2] = sp.absolute(sp.average(f1[:, 2])-sp.average(f2[:, 2]))
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
            logger.warning('The supplied pores are not coplanar. Area will be'
                           'approximate')
            pass
        return A

    def _compress_labels(self, label_array):
        # Make cluster number contiguous
        array = sp.array(label_array)
        if array.dtype != int:
            raise Exception('label_array must be intergers')
        min_val = sp.amin(array)
        if min_val >= 0:
            min_val = 0
        array = array + sp.absolute(min_val)
        nums = sp.unique(array)
        temp = sp.zeros((sp.amax(array)+1,))
        temp[nums] = sp.arange(0, sp.size(nums))
        array = temp[array].astype(array.dtype)
        return array

    def add_boundary_pores(self, pores, offset, apply_label='boundary'):
        r"""
        This method uses ``clone_pores`` to clone the input pores, then shifts
        them the specified amount and direction, then applies the given label.

        Parameters
        ----------
        pores : array_like
            List of pores to offset.  If no pores are specified, then it
            assumes that all surface pores are to be cloned.

        offset : 3 x 1 array
            The distance in vector form which the cloned boundary pores should
            be offset.  If no spacing is provided, then the spacing is inferred
            from the Network.

        apply_label : string
            This label is applied to the boundary pores.  Default is
            'boundary'.

        Examples
        --------
        >>> import OpenPNM as op
        >>> pn = op.Network.Cubic(shape=[5, 5, 5])
        >>> print(pn.Np)  # Confirm initial Network size
        125
        >>> Ps = pn.pores('top')  # Select pores on top face
        >>> pn.add_boundary_pores(pores=Ps, offset=[0, 0, 1])
        >>> print(pn.Np)  # Confirm addition of 25 new pores
        150
        >>> 'pore.boundary' in pn.labels()  # Default label is created
        True
        """
        # Parse the input pores
        Ps = sp.array(pores, ndmin=1)
        if Ps.dtype is bool:
            Ps = self.toindices(Ps)
        if sp.size(pores) == 0:  # Handle an empty array if given
            return sp.array([], dtype=sp.int64)
        # Clone the specifed pores
        self.clone_pores(pores=Ps)
        newPs = self.pores('pore.clone')
        del self['pore.clone']
        # Offset the cloned pores
        self['pore.coords'][newPs] += offset
        # Apply labels to boundary pores (trim leading 'pores' if present)
        label = apply_label.split('.')[-1]
        label = 'pore.' + label
        logger.debug('The label \''+label+'\' has been applied')
        self[label] = False
        self[label][newPs] = True
