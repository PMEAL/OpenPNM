import itertools
import uuid
import scipy as sp
import scipy.sparse as sprs
import scipy.spatial as sptl
from openpnm.core import Base, Simulation, Workspace, ModelsMixin, logging
logger = logging.getLogger()
ws = Workspace()


class GenericNetwork(Base, ModelsMixin):
    r"""
    GenericNetwork - Base class to construct pore networks

    Parameters
    ----------
    name : string
        Unique name for Network object

    """

    _prefix = 'net'

    def __init__(self, simulation=None, **kwargs):
        if simulation is None:
            simulation = Simulation()
        super().__init__(simulation=simulation, **kwargs)
        self['pore._id'] = [str(uuid.uuid4()) for i in self.Ps]
        self['throat._id'] = [str(uuid.uuid4()) for i in self.Ts]

        # Initialize adjacency and incidence matrix dictionaries
        self._im = {}
        self._am = {}

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
        super().__setitem__(prop, value)

    def __getitem__(self, key):
        if key.split('.')[-1] == self.name:
            element = key.split('.')[0]
            return self[element+'.all']
        if key not in self.keys():
            logger.debug(key + ' not on Network, constructing data from Geometries')
            return self._interleave_data(key, self.simulation.geometries.values())
        else:
            return super().__getitem__(key)

    def get_adjacency_matrix(self, fmt='coo'):
        # Retrieve existing matrix if available
        if fmt in self._am.keys():
            return self._am[fmt]
        else:
            am = self.create_adjacency_matrix(fmt=fmt)
        return am

    def get_incidence_matrix(self, fmt='coo'):
        if fmt in self._im.keys():
            return self._im[fmt]
        im = self.create_incidence_matrix(fmt=fmt)
        return im

    im = property(fget=get_incidence_matrix)

    am = property(fget=get_adjacency_matrix)

    def create_adjacency_matrix(self, data=None, fmt='coo'):
        r"""
        Generates a weighted adjacency matrix in the desired sparse format

        Parameters
        ----------
        data : array_like, optional
            An Nt-long array containing the throat values to enter into the
            matrix (in graph theory these are known as the 'weights').  If
            omitted, ones are used to create a standard adjacency matrix
            representing connectivity only.  Zero values are dropped from the
            matrix.

        fmt : string, optional
            The sparse storage format to return.  Options are:

            **'coo'** : (default) This is the native format of OpenPNM data

            **'lil'** : Enables row-wise slice of data

            **'csr'** : Favored by most linear algebra routines

        Returns
        -------
        Returns an adjacency matrix in the specified Scipy sparse format

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> vals = sp.rand(pn.num_throats(),) < 0.5
        >>> temp = pn.create_adjacency_matrix(data=vals, fmt='csr')

        """
        # Check if provided data is valid
        data_flag = False
        if data is None:
            data = sp.ones((self.Nt,))
            data_flag = True
        elif sp.shape(data)[0] != self.Nt:
            raise Exception('Received weights are incorrect length')

        # Append row & col to each other, and data to itself
        ind = data != 0
        conn = self['throat.conns'][ind]
        row = conn[:, 0]
        col = conn[:, 1]
        data = data[ind]
        row = sp.append(row, conn[:, 1])
        col = sp.append(col, conn[:, 0])
        data = sp.append(data, data)

        # Generate sparse adjacency matrix in 'coo' format
        temp = sprs.coo_matrix((data, (row, col)), (self.Np, self.Np))

        # Save the 'clean' matrix on object for future use
        if data_flag:
            if 'coo' in self._am.keys():
                self._am['coo'] = temp
            if 'csr' in self._am.keys():
                self._am['csr'] = temp.tocsr()
            if 'lil' in self._am.keys():
                self._am['lil'] = temp.tolil()

        # Convert to requested format
        if fmt == 'coo':
            pass  # temp is already in coo format
        if fmt == 'csr':
            temp = temp.tocsr()
        if fmt == 'lil':
            temp = temp.tolil()

        return temp

    def create_incidence_matrix(self, data=None, fmt='coo'):
        r"""
        Creates an incidence matrix filled with supplied throat values

        Parameters
        ----------
        data : array_like, optional
            An array containing the throat values to enter into the matrix (In
            graph theory these are known as the 'weights').  If omitted, ones
            are used to create a standard incidence matrix representing
            connectivity only.  Zero values are dropped from the matrix.

        fmt : string, optional
            The sparse storage format to return.  Options are:

            **'coo'** : (default) This is the native format of OpenPNMs data

            **'lil'** : Enables row-wise slice of data

            **'csr'** : Favored by most linear algebra routines

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
        # Check if provided data is valid
        data_flag = False
        if data is None:
            data = sp.ones((self.Nt,))
            data_flag
        elif sp.shape(data)[0] != self.Nt:
            raise Exception('Received dataset of incorrect length')

        ind = data > 0
        conn = self['throat.conns'][ind]
        row = conn[:, 0]
        row = sp.append(row, conn[:, 1])
        col = self.throats('all')[ind]
        col = sp.append(col, col)
        data = sp.append(data[ind], data[ind])

        temp = sprs.coo.coo_matrix((data, (row, col)), (self.Np, self.Nt))

        # Save the 'clean' matrix on object for future use
        if data_flag:
            if 'coo' in self._am.keys():
                self._im['coo'] = temp
            if 'csr' in self._am.keys():
                self._im['csr'] = temp.tocsr()
            if 'lil' in self._am.keys():
                self._im['lil'] = temp.tolil()

        # Convert to requested format
        if fmt == 'coo':
            pass  # temp is already in coo format
        if fmt == 'csr':
            temp = temp.tocsr()
        if fmt == 'lil':
            temp = temp.tolil()

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
        1D array (if flatten is True) or ndarray of arrays (if flatten is
        False)

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
        Ts = self._parse_indices(throats)
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
        in the Numpy sense, so could be slow with large P1, P2 inputs
        """
        P1 = self._parse_indices(P1)
        P2 = self._parse_indices(P2)
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
        ``flatten`` is ``False``)

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
        pores = self._parse_indices(pores)
        if sp.size(pores) == 0:
            return sp.array([], ndmin=1, dtype=int)

        # Test for existence of incidence or adjacency matrix
        if element == 'pore':
            temp = self.get_adjacency_matrix(fmt='lil')
            neighbors = temp.rows[[pores]]
        elif element == 'throat':
            temp = self.get_incidence_matrix(fmt='lil')
            neighbors = temp.rows[[pores]]

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
            neighbors = sp.array(neighbors, ndmin=1, dtype=int)
        else:
            # Convert lists in array to numpy arrays
            neighbors = [sp.array(neighbors[i]) for i in range(0, len(pores))]
            neighbors = sp.array(neighbors, ndmin=1)
        return neighbors

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
        pores = self._parse_indices(pores)
        # Count number of neighbors
        num = self._find_neighbors(pores, element=element, flatten=flatten,
                                   mode=mode, excl_self=True)
        num = sp.array([sp.size(i) for i in num], dtype=int)
        if flatten:
            num = sp.sum(num)
            num = int(num)
        return num

    def find_nearby_pores(self, pores, R, flatten=False, excl_self=True):
        r"""
        Find all pores within a given radial distance of the input pore(s)
        regardless of whether or not they are toplogically connected.

        Parameters
        ----------
        pores : array_like
            The list of pores for whom nearby neighbors are to be found

        R : scalar
            The maximum radius within which the search should be performed

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
        pores = self._parse_indices(pores)
        # Handle an empty array if given
        if sp.size(pores) == 0:
            return sp.array([], dtype=sp.int64)
        if R <= 0:
            logger.error('Provided distances should be greater than 0')
            if flatten:
                Pn = sp.array([])
            else:
                Pn = sp.array([sp.array([]) for i in range(0, len(pores))])
            return Pn.astype(sp.int64)
        # Create kdTree objects
        kd = sptl.cKDTree(self['pore.coords'])
        kd_pores = sptl.cKDTree(self['pore.coords'][pores])
        # Perform search
        Pn = kd_pores.query_ball_tree(kd, r=R)
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
