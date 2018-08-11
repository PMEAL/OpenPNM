import uuid
import scipy as sp
import scipy.sparse as sprs
import scipy.spatial as sptl
import scipy.sparse.csgraph as csg
from openpnm.core import Base, ModelsMixin
from openpnm import topotools
from openpnm.utils import HealthDict, Workspace, logging
logger = logging.getLogger()
ws = Workspace()


class GenericNetwork(Base, ModelsMixin):
    r"""
    This generic class contains the main functionality used by all networks

    Parameters
    ----------
    coords : array_like
        An Np-by-3 array of [x, y, z] coordinates for each pore.

    conns : array_like
        An Nt-by-2 array of [head, tail] connections between pores.

    Notes
    -----
    The GenericNetwork class houses a number of methods used for querying and
    managing the network's spatial and topological information.  The following
    table gives a very short overview of the methods added those already found
    on the ``openpnm.core.Base`` class.

    +-----------------------------+-------------------------------------------+
    | Method or Attribute         | Functionality                             |
    +=============================+===========================================+
    | ``create_adjacency_matrix`` | Create an adjacency matrix using given    |
    |                             | weights in a specified format             |
    +-----------------------------+-------------------------------------------+
    | ``create_incidence_matrix`` | Create an incidence matrix using given    |
    |                             | weights in a specified format             |
    +-----------------------------+-------------------------------------------+
    | ``get_adjacency_matrix``    | Retrieve an existing adjacency matrix in  |
    |                             | the specified format (from ``am``)        |
    +-----------------------------+-------------------------------------------+
    | ``get_incidence_matrix``    | Retrieve an existing incidence matrix in  |
    |                             | the specified format (from ``im``)        |
    +-----------------------------+-------------------------------------------+
    | ``am``                      | Returns the adjacency matrix in COO format|
    +-----------------------------+-------------------------------------------+
    | ``im``                      | Returns the incidence matrix in COO format|
    +-----------------------------+-------------------------------------------+
    | ``find_neighbor_pores``     | For a given set of pores, find all        |
    |                             | neighboring pores                         |
    +-----------------------------+-------------------------------------------+
    | ``find_neighbor_throats``   | For a given set of pores, find all        |
    |                             | neighboring throats                       |
    +-----------------------------+-------------------------------------------+
    | ``find_connecting_throat``  | For each pair of throats find the pores   |
    |                             | they connect                              |
    +-----------------------------+-------------------------------------------+
    | ``find_connected_pores``    | For each throat, find the pores which it  |
    |                             | connects                                  |
    +-----------------------------+-------------------------------------------+
    | ``num_neighbors``           | For a given set of pores find the number  |
    |                             | of neighbors for each                     |
    +-----------------------------+-------------------------------------------+
    | ``find_nearby_pores``       | For a given set of pores, find pores that |
    |                             | are within a certain distance             |
    +-----------------------------+-------------------------------------------+
    | ``check_network_health``    | Check the topology for any problems such  |
    |                             | as isolated pores                         |
    +-----------------------------+-------------------------------------------+

    Examples
    --------
    >>> import openpnm as op

    Create some pore coordinates and connections manually and assign to a
    GenericNetwork instance.  Consider a linear network of 4 pores and 3
    throats:

    ::

        0 ―― 1 ―― 3 ―― 2

    >>> coords = [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]]
    >>> conns = [[0, 1], [1, 3], [2, 3]]
    >>> pn = op.network.GenericNetwork(conns=conns, coords=coords)

    Networks have two required properties: 'pore.coords' and 'throat.conns'.
    These arrays indicate the spatial location of each pore, and which pores
    are connected to which.  Without these the Network object cannot function.

    >>> print(pn.props())
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    1     : pore.coords
    2     : throat.conns
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――

    The GenericNetwork class has several methods for querying the topology.

    >>> Ps = pn.find_neighbor_pores(pores=1)
    >>> print(Ps)
    [0 3]
    >>> Ts = pn.find_neighbor_throats(pores=[0, 1])
    >>> print(Ts)
    [0 1]
    >>> print(pn.num_neighbors(2))
    [1]

    All of the topological queries are accomplished by inspecting the adjacency
    and incidence matrices.  They are created on demand, and are stored for
    future use to save construction time.

    """
    def __init__(self, conns=None, coords=None, project=None, settings={},
                 **kwargs):
        self.settings.setdefault('prefix', 'net')
        self.settings.update(settings)
        super().__init__(project=project, **kwargs)
        if coords is not None:
            Np = sp.shape(coords)[0]
            self['pore.all'] = sp.ones(Np, dtype=bool)
            self['pore.coords'] = sp.array(coords)
        if conns is not None:
            Nt = sp.shape(conns)[0]
            self['throat.all'] = sp.ones(Nt, dtype=bool)
            self['throat.conns'] = sp.array(conns)
        # Initialize adjacency and incidence matrix dictionaries
        self._im = {}
        self._am = {}

    def __setitem__(self, key, value):
        if key == 'throat.conns':
            if sp.shape(value)[1] != 2:
                logger.error('Wrong size for throat conns!')
            else:
                if sp.any(value[:, 0] > value[:, 1]):
                    logger.warning('Converting throat.conns to be upper ' +
                                   'triangular')
                    value = sp.sort(value, axis=1)
        if self.project:
            for item in self.project.geometries().values():
                exclude = {'pore.all', 'throat.all'}
                if key in set(item.keys()).difference(exclude):
                    raise Exception(key+' already exists on '+item.name)
        super().__setitem__(key, value)

    def __getitem__(self, key):
        # Deal with special keys first
        if key.split('.')[-1] == self.name:
            element = key.split('.')[0]
            return self[element+'.all']
        if key.split('.')[-1] == '_id':
            self._gen_ids()
        # Now get values if present, or regenerate them
        vals = self.get(key)
        if vals is None:  # Invoke interleave data
            vals = self.interleave_data(key)
        return vals

    def _gen_ids(self):
        IDs = super().get('pore._id', sp.array([], ndmin=1, dtype=sp.int64))
        if len(IDs) < self.Np:
            temp = ws._gen_ids(size=self.Np - len(IDs))
            self['pore._id'] = sp.concatenate((IDs, temp))
        IDs = super().get('throat._id', sp.array([], ndmin=1, dtype=sp.int64))
        if len(IDs) < self.Nt:
            temp = ws._gen_ids(size=self.Nt - len(IDs))
            self['throat._id'] = sp.concatenate((IDs, temp))

    def get_adjacency_matrix(self, fmt='coo'):
        r"""
        Returns an adjacency matrix in the specified sparse format, with 1's
        indicating the non-zero values.

        Parameters
        ----------
        fmt : string, optional
            The sparse storage format to return.  Options are:

            **'coo'** : (default) This is the native format of OpenPNM data

            **'lil'** : Enables row-wise slice of the matrix

            **'csr'** : Favored by most linear algebra routines

            **'dok'** : Enables subscript access of locations

        Notes
        -----
        This method will only create the requested matrix in the specified
        format if one is not already saved on the object.  If not present,
        this method will create and return the matrix, as well as store it
        for future use.

        To obtain a matrix with weights other than ones at each non-zero
        location use ``create_adjacency_matrix``.
        """
        # Retrieve existing matrix if available
        if fmt in self._am.keys():
            am = self._am[fmt]
        elif self._am.keys():
            am = self._am[list(self._am.keys())[0]]
            tofmt = getattr(am, 'to'+fmt)
            am = tofmt()
            self._am[fmt] = am
        else:
            am = self.create_adjacency_matrix(weights=self.Ts, fmt=fmt)
            self._am[fmt] = am
        return am

    def get_incidence_matrix(self, fmt='coo'):
        r"""
        Returns an incidence matrix in the specified sparse format, with 1's
        indicating the non-zero values.

        Parameters
        ----------
        fmt : string, optional
            The sparse storage format to return.  Options are:

            **'coo'** : (default) This is the native format of OpenPNM data

            **'lil'** : Enables row-wise slice of the matrix

            **'csr'** : Favored by most linear algebra routines

            **'dok'** : Enables subscript access of locations

        Notes
        -----
        This method will only create the requested matrix in the specified
        format if one is not already saved on the object.  If not present,
        this method will create and return the matrix, as well as store it
        for future use.

        To obtain a matrix with weights other than ones at each non-zero
        location use ``create_incidence_matrix``.
        """
        if fmt in self._im.keys():
            im = self._im[fmt]
        elif self._im.keys():
            im = self._am[list(self._im.keys())[0]]
            tofmt = getattr(im, 'to'+fmt)
            im = tofmt()
            self._im[fmt] = im
        else:
            im = self.create_incidence_matrix(weights=self.Ts, fmt=fmt)
            self._im[fmt] = im
        return im

    im = property(fget=get_incidence_matrix)

    am = property(fget=get_adjacency_matrix)

    def create_adjacency_matrix(self, weights=None, fmt='coo', triu=False,
                                drop_zeros=False):
        r"""
        Generates a weighted adjacency matrix in the desired sparse format

        Parameters
        ----------
        weights : array_like, optional
            An array containing the throat values to enter into the matrix
            (in graph theory these are known as the 'weights').

            If the array is Nt-long, it implies that the matrix is symmetric,
            so the upper and lower triangular regions are mirror images.  If
            it is 2*Nt-long then it is assumed that the first Nt elements are
            for the upper triangle, and the last Nt element are for the lower
            triangular.

            If omitted, ones are used to create a standard adjacency matrix
            representing connectivity only.

        fmt : string, optional
            The sparse storage format to return.  Options are:

            **'coo'** : (default) This is the native format of OpenPNM data

            **'lil'** : Enables row-wise slice of the matrix

            **'csr'** : Favored by most linear algebra routines

            **'dok'** : Enables subscript access of locations

        triu : boolean (default is ``False``)
            If ``True``, the returned sparse matrix only contains the upper-
            triangular elements.  This argument is ignored if the ``weights``
            array is 2*Nt-long.

        drop_zeros : boolean (default is ``False``)
            If ``True``, applies the ``eliminate_zeros`` method of the sparse
            array to remove all zero locations.

        Returns
        -------
        An adjacency matrix in the specified Scipy sparse format.

        Notes
        -----
        The adjacency matrix is used by OpenPNM for finding the pores
        connected to a give pore or set of pores.  Specifically, an adjacency
        matrix has Np rows and Np columns.  Each row represents a pore,
        containing non-zero values at the locations corresponding to the
        indices of the pores connected to that pore.  The ``weights`` argument
        indicates what value to place at each location, with the default
        being 1's to simply indicate connections. Another useful option is
        throat indices, such that the data values on each row indicate which
        throats are connected to the pore.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> weights = sp.rand(pn.num_throats(), ) < 0.5
        >>> am = pn.create_adjacency_matrix(weights=weights, fmt='csr')

        """
        # Check if provided data is valid
        if weights is None:
            weights = sp.ones((self.Nt,), dtype=int)
        elif sp.shape(weights)[0] not in [self.Nt, 2*self.Nt, (self.Nt, 2)]:
            raise Exception('Received weights are of incorrect length')

        # Append row & col to each other, and data to itself
        conn = self['throat.conns']
        row = conn[:, 0]
        col = conn[:, 1]
        if weights.shape == (2*self.Nt,):
            row = sp.append(row, conn[:, 1])
            col = sp.append(col, conn[:, 0])
        elif weights.shape == (self.Nt, 2):
            row = sp.append(row, conn[:, 1])
            col = sp.append(col, conn[:, 0])
            weights = weights.flatten(order='F')
        elif not triu:
            row = sp.append(row, conn[:, 1])
            col = sp.append(col, conn[:, 0])
            weights = sp.append(weights, weights)

        # Generate sparse adjacency matrix in 'coo' format
        temp = sprs.coo_matrix((weights, (row, col)), (self.Np, self.Np))

        if drop_zeros:
            temp.eliminate_zeros()

        # Convert to requested format
        if fmt == 'coo':
            pass  # temp is already in coo format
        elif fmt == 'csr':
            temp = temp.tocsr()
        elif fmt == 'lil':
            temp = temp.tolil()
        elif fmt == 'dok':
            temp = temp.todok()

        return temp

    def create_incidence_matrix(self, weights=None, fmt='coo',
                                drop_zeros=False):
        r"""
        Creates a weighted incidence matrix in the desired sparse format

        Parameters
        ----------
        weights : array_like, optional
            An array containing the throat values to enter into the matrix (In
            graph theory these are known as the 'weights').  If omitted, ones
            are used to create a standard incidence matrix representing
            connectivity only.

        fmt : string, optional
            The sparse storage format to return.  Options are:

            **'coo'** : (default) This is the native format of OpenPNMs data

            **'lil'** : Enables row-wise slice of the matrix

            **'csr'** : Favored by most linear algebra routines

            **'dok'** : Enables subscript access of locations

        drop_zeros : boolean (default is ``False``)
            If ``True``, applies the ``eliminate_zeros`` method of the sparse
            array to remove all zero locations.

        Returns
        -------
        An incidence matrix in the specified sparse format

        Notes
        -----
        The incidence matrix is a cousin to the adjacency matrix, and used by
        OpenPNM for finding the throats connected to a give pore or set of
        pores.  Specifically, an incidence matrix has Np rows and Nt columns,
        and each row represents a pore, containing non-zero values at the
        locations corresponding to the indices of the throats connected to that
        pore.  The ``weights`` argument indicates what value to place at each
        location, with the default being 1's to simply indicate connections.
        Another useful option is throat indices, such that the data values
        on each row indicate which throats are connected to the pore, though
        this is redundant as it is identical to the locations of non-zeros.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> weights = sp.rand(pn.num_throats(), ) < 0.5
        >>> im = pn.create_incidence_matrix(weights=weights, fmt='csr')
        """
        # Check if provided data is valid
        if weights is None:
            weights = sp.ones((self.Nt,), dtype=int)
        elif sp.shape(weights)[0] != self.Nt:
            raise Exception('Received dataset of incorrect length')

        conn = self['throat.conns']
        row = conn[:, 0]
        row = sp.append(row, conn[:, 1])
        col = sp.arange(self.Nt)
        col = sp.append(col, col)
        weights = sp.append(weights, weights)

        temp = sprs.coo.coo_matrix((weights, (row, col)), (self.Np, self.Nt))

        if drop_zeros:
            temp.eliminate_zeros()

        # Convert to requested format
        if fmt == 'coo':
            pass  # temp is already in coo format
        elif fmt == 'csr':
            temp = temp.tocsr()
        elif fmt == 'lil':
            temp = temp.tolil()
        elif fmt == 'dok':
            temp = temp.todok()

        return temp

    def find_connected_pores(self, throats=[], flatten=False, mode='union'):
        r"""
        Return a list of pores connected to the given list of throats

        Parameters
        ----------
        throats : array_like
            List of throats numbers

        flatten : boolean, optional
            If ``True`` (default) a 1D array of unique pore numbers is
            returned. If ``False`` each location in the the returned array
            contains a sub-arras of neighboring pores for each input throat,
            in the order they were sent.

        mode : string
            Specifies logic to filter the resulting list.  Options are:

            **'or'** : (default) All neighbors of the input throats.  This is
            also known as the 'union' in set theory or 'any' in boolean logic.
            Both keywords are accepted and treated as 'or'.

            **'xor'** : Only neighbors of one and only one input throat.  This
            is useful for finding the sites that are not shared by any of the
            input throats.

            **'xnor'** : Neighbors that are shared by two or more input
            throats. This is equivalent to finding all neighbors with 'or',
            minus those found with 'xor', and is useful for finding neighbors
            that the inputs have in common.

            **'and'** : Only neighbors shared by all input throats.  This is
            also known as 'intersection' in set theory and (somtimes) as 'all'
            in boolean logic.  Both keywords are accepted and treated as 'and'.

        Returns
        -------
        1D array (if ``flatten`` is ``True``) or ndarray of arrays (if
        ``flatten`` is ``False``)

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.find_connected_pores(throats=[0, 1])
        array([[0, 1],
               [1, 2]])
        >>> pn.find_connected_pores(throats=[0, 1], flatten=True)
        array([0, 1, 2])

        """
        Ts = self._parse_indices(throats)
        am = self.get_adjacency_matrix(fmt='coo')
        pores = topotools.find_connected_sites(bonds=Ts, am=am,
                                               flatten=flatten, logic=mode)
        return pores

    def find_connecting_throat(self, P1, P2):
        r"""
        Return the throat index connecting pairs of pores

        Parameters
        ----------
        P1 , P2 : array_like
            The indices of the pores whose throats are sought.  These can be
            vectors of indices, but must be the same length

        Returns
        -------
        Returns a list the same length as P1 (and P2) with the each element
        containing the throat index that connects the corresponding pores,
        or `None`` if pores are not connected.

        Notes
        -----
        The returned list can be converted to an ND-array, which will convert
        the ``None`` values to ``nan``.  These can then be found using
        ``scipy.isnan``.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.find_connecting_throat([0, 1, 2], [2, 2, 2])
        [None, 1, None]
        """
        am = self.create_adjacency_matrix(weights=self.Ts, fmt='coo')
        sites = sp.vstack((P1, P2)).T
        Ts = topotools.find_connecting_bonds(sites=sites, am=am)
        return Ts

    def find_neighbor_pores(self, pores, mode='union', flatten=True,
                            include_input=False):
        r"""
        Returns a list of pores that are direct neighbors to the given pore(s)

        Parameters
        ----------
        pores : array_like
            Indices of the pores whose neighbors are sought

        flatten : boolean
            If ``True`` (default) the returned result is a compressed array of
            all neighbors.  If ``False``, a list of lists with each sub-list
            containing the neighbors for each input site.  Note that an
            *unflattened* list might be slow to generate since it is a Python
            ``list`` rather than a Numpy ``array``.

        include_input : bool
            If ``False`` (default) then the input pores are not included in
            the returned list(s).

        mode : string
            Specifies logic to filter the resulting list.  Options are:

            **'or'** : (default) All neighbors of the input pores.  This is
            also known as the 'union' in set theory or 'any' in boolean logic.
            Both keywords are accepted and treated as 'or'.

            **'xor'** : Only neighbors of one and only one input pore.  This
            is useful for finding the pores that are not shared by any of the
            input pores.  This is known as 'exclusive_or' in set theory, and
            is an accepted input.

            **'xnor'** : Neighbors that are shared by two or more input pores.
            This is equivalent to finding all neighbors with 'or', minus those
            found with 'xor', and is useful for finding neighbors that the
            inputs have in common.

            **'and'** : Only neighbors shared by all input pores.  This is also
            known as 'intersection' in set theory and (somtimes) as 'all' in
            boolean logic.  Both keywords are accepted and treated as 'and'.

        Returns
        -------
        If ``flatten`` is ``True``, returns a 1D array of pore indices filtered
        according to the specified mode.  If ``flatten`` is ``False``, returns
        a list of lists, where each list contains the neighbors of the
        corresponding input pores.

        Notes
        -----
        The ``logic`` options are applied to neighboring pores only, thus it
        is not possible to obtain pores that are part of the global set but
        not neighbors. This is because (a) the list of global pores might be
        very large, and (b) it is not possible to return a list of neighbors
        for each input pores if global pores are considered.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.find_neighbor_pores(pores=[0, 2])
        array([ 1,  3,  5,  7, 25, 27])
        >>> pn.find_neighbor_pores(pores=[0, 1])
        array([ 2,  5,  6, 25, 26])
        >>> pn.find_neighbor_pores(pores=[0, 1], mode='union',
        ...                        include_input=True)
        array([ 0,  1,  2,  5,  6, 25, 26])
        >>> pn.find_neighbor_pores(pores=[0, 2], flatten=False)
        [[1, 5, 25], [1, 3, 7, 27]]
        >>> pn.find_neighbor_pores(pores=[0, 2], mode='intersection')
        array([1])
        >>> pn.find_neighbor_pores(pores=[0, 2], mode='exclusive_or')
        array([ 3,  5,  7, 25, 27])
        """
        pores = self._parse_indices(pores)
        if sp.size(pores) == 0:
            return sp.array([], ndmin=1, dtype=int)
        if 'lil' not in self._am.keys():
            self.get_adjacency_matrix(fmt='lil')
        neighbors = topotools.find_neighbor_sites(sites=pores, logic=mode,
                                                  am=self._am['lil'],
                                                  flatten=flatten,
                                                  include_input=include_input)
        return neighbors

    def find_neighbor_throats(self, pores, mode='union', flatten=True):
        r"""
        Returns a list of throats neighboring the given pore(s)

        Parameters
        ----------
        pores : array_like
            Indices of pores whose neighbors are sought

        flatten : boolean, optional
            If ``True`` (default) a 1D array of unique throat indices is
            returned. If ``False`` the returned array contains arrays of
            neighboring throat indices for each input pore, in the order
            they were sent.

        mode : string
            Specifies logic to filter the resulting list.  Options are:

            **'or'** : (default) All neighbors of the input pores.  This is
            also known as the 'union' in set theory or 'any' in boolean logic.
            Both keywords are accepted and treated as 'or'.

            **'xor'** : Only neighbors of one and only one input pore.  This
            is useful for finding the thraots that are not shared by any of the
            input pores.

            **'xnor'** : Neighbors that are shared by two or more input pores.
            This is equivalent to finding all neighbors with 'or', minus those
            found with 'xor', and is useful for finding neighbors that the
            inputs have in common.

            **'and'** : Only neighbors shared by all input pores.  This is also
            known as 'intersection' in set theory and (somtimes) as 'all' in
            boolean logic.  Both keywords are accepted and treated as 'and'.

        Returns
        -------
        If ``flatten`` is ``True``, returns a 1D array of throat indices
        filtered according to the specified mode.  If ``flatten`` is ``False``,
        returns a list of lists, where each list contains the neighbors of the
        corresponding input pores.

        Notes
        -----
        The ``logic`` options are applied to neighboring bonds only, thus it
        is not possible to obtain bonds that are part of the global set but
        not neighbors. This is because (a) the list of global bonds might be
        very large, and (b) it is not possible to return a list of neighbors
        for each input site if global sites are considered.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.find_neighbor_throats(pores=[0, 1])
        array([  0,   1, 100, 101, 200, 201])
        >>> pn.find_neighbor_throats(pores=[0, 1], flatten=False)
        [[0, 100, 200], [0, 1, 101, 201]]

        """
        pores = self._parse_indices(pores)
        if sp.size(pores) == 0:
            return sp.array([], ndmin=1, dtype=int)
        if 'lil' not in self._im.keys():
            self.get_incidence_matrix(fmt='lil')
        neighbors = topotools.find_neighbor_bonds(sites=pores, logic=mode,
                                                  im=self._im['lil'],
                                                  flatten=flatten)
        return neighbors

    def _find_neighbors(self, pores, element, **kwargs):
        element = self._parse_element(element=element, single=True)
        if sp.size(pores) == 0:
            return sp.array([], ndmin=1, dtype=int)
        if element == 'pore':
            neighbors = self.find_neighbor_pores(pores=pores, **kwargs)
        else:
            neighbors = self.find_neighbor_throats(pores=pores, **kwargs)
        return neighbors

    def num_neighbors(self, pores, mode='or', flatten=False):
        r"""
        Returns the number of neigbhoring pores for each given input pore

        Parameters
        ----------
        pores : array_like
            Pores whose neighbors are to be counted

        flatten : boolean (optional)
            If ``False`` (default) the number of pores neighboring each input
            pore as an array the same length as ``pores``.  If ``True`` the
            sum total number of is counted.

        mode : string
            The logic to apply to the returned count of pores.

            **'or'** : (default) All neighbors of the input pores.  This is
            also known as the 'union' in set theory or 'any' in boolean logic.
            Both keywords are accepted and treated as 'or'.

            **'xor'** : Only neighbors of one and only one input pore.  This
            is useful for counting the pores that are not shared by any of the
            input pores.  This is known as 'exclusive_or' in set theory, and
            is an accepted input.

            **'xnor'** : Neighbors that are shared by two or more input pores.
            This is equivalent to counting all neighbors with 'or', minus those
            found with 'xor', and is useful for finding neighbors that the
            inputs have in common.

            **'and'** : Only neighbors shared by all input pores.  This is also
            known as 'intersection' in set theory and (somtimes) as 'all' in
            boolean logic.  Both keywords are accepted and treated as 'and'.

        Returns
        -------
        If ``flatten`` is False, a 1D array with number of neighbors in each
        element, otherwise a scalar value of the number of neighbors.

        Notes
        -----
        This method literally just counts the number of elements in the array
        returned by ``find_neighbor_pores`` using the same logic.  Explore
        those methods if uncertain about the meaning of the ``mode`` argument
        here.

        See Also
        --------
        find_neighbor_pores
        find_neighbor_throats

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.num_neighbors(pores=[0, 1], flatten=False)
        array([3, 4])
        >>> pn.num_neighbors(pores=[0, 2], flatten=True)
        6
        >>> pn.num_neighbors(pores=[0, 2], mode='and', flatten=True)
        1
        """
        pores = self._parse_indices(pores)
        # Count number of neighbors
        num = self.find_neighbor_pores(pores, flatten=flatten,
                                       mode=mode, include_input=False)
        if flatten:
            num = sp.size(num)
        else:
            num = sp.array([sp.size(i) for i in num], dtype=int)
        return num

    def find_nearby_pores(self, pores, r, flatten=False, include_input=False):
        r"""
        Find all pores within a given radial distance of the input pore(s)
        regardless of whether or not they are toplogically connected.

        Parameters
        ----------
        pores : array_like
            The list of pores for whom nearby neighbors are to be found

        r : scalar
            The maximum radius within which the search should be performed

        include_input : bool
            Controls whether the input pores should be included in the returned
            list.  The default is ``False`` which means they are not included.

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
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[3, 3, 3])
        >>> pn.find_nearby_pores(pores=[0, 1], r=1)
        array([array([1, 3, 9]), array([ 0,  2,  4, 10])], dtype=object)
        >>> pn.find_nearby_pores(pores=[0, 1], r=0.5)
        array([], shape=(2, 0), dtype=int64)
        >>> pn.find_nearby_pores(pores=[0, 1], r=1, flatten=True)
        array([ 2,  3,  4,  9, 10])
        """
        pores = self._parse_indices(pores)
        # Handle an empty array if given
        if sp.size(pores) == 0:
            return sp.array([], dtype=sp.int64)
        if r <= 0:
            raise Exception('Provided distances should be greater than 0')
        # Create kdTree objects
        kd = sptl.cKDTree(self['pore.coords'])
        kd_pores = sptl.cKDTree(self['pore.coords'][pores])
        # Perform search
        Pn = kd_pores.query_ball_tree(kd, r=r)
        # Sort the indices in each list
        # [Pn[i].sort() for i in range(0, sp.size(pores))]
        if flatten:  # Convert list of lists to a flat nd-array
            temp = sp.concatenate((Pn))
            Pn = sp.unique(temp)
            if include_input:  # Remove inputs if necessary
                Pn = Pn[~sp.in1d(Pn, pores)]
        else:  # Convert list of lists to an nd-array of nd-arrays
            if include_input:  # Remove inputs if necessary
                [Pn[i].remove(pores[i]) for i in range(0, sp.size(pores))]
            temp = [sp.array(Pn[i]) for i in range(0, sp.size(pores))]
            Pn = sp.array(temp)
        if Pn.dtype == float:
            Pn = Pn.astype(sp.int64)
        return Pn

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
        - This is just a 'check' and does not 'fix' the problems it finds
        """

        health = HealthDict()
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
            return health

        # Check for throats that loop back onto the same pore
        P12 = self['throat.conns']
        hits = sp.where(P12[:, 0] == P12[:, 1])[0]
        if sp.size(hits) > 0:
            health['looped_throats'] = hits

        # Check for individual isolated pores
        Ps = self.num_neighbors(self.pores())
        if sp.sum(Ps == 0) > 0:
            health['isolated_pores'] = sp.where(Ps == 0)[0]

        # Check for separated clusters of pores
        temp = []
        am = self.create_adjacency_matrix(fmt='coo', triu=True)
        Cs = csg.connected_components(am, directed=False)[1]
        if sp.unique(Cs).size > 1:
            for i in sp.unique(Cs):
                temp.append(sp.where(Cs == i)[0])
            b = sp.array([len(item) for item in temp])
            c = sp.argsort(b)[::-1]
            for i in range(0, len(c)):
                health['disconnected_clusters'].append(temp[c[i]])
                if i > 0:
                    health['trim_pores'].extend(temp[c[i]])

        # Check for duplicate throats
        am = self.create_adjacency_matrix(fmt='csr', triu=True).tocoo()
        hits = sp.where(am.data > 1)[0]
        if len(hits):
            mergeTs = []
            hits = sp.vstack((am.row[hits], am.col[hits])).T
            ihits = hits[:, 0] + 1j*hits[:, 1]
            conns = self['throat.conns']
            iconns = conns[:, 0] + 1j*conns[:, 1]  # Convert to imaginary
            for item in ihits:
                mergeTs.append(sp.where(iconns == item)[0])
            health['duplicate_throats'] = mergeTs

        # Check for bidirectional throats
        adjmat = self.create_adjacency_matrix(fmt='coo')
        num_full = adjmat.sum()
        temp = sprs.triu(adjmat, k=1)
        num_upper = temp.sum()
        if num_full > num_upper:
            biTs = sp.where(self['throat.conns'][:, 0] >
                            self['throat.conns'][:, 1])[0]
            health['bidirectional_throats'] = biTs.tolist()

        return health
