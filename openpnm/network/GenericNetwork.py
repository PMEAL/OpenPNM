import uuid
import scipy as sp
import scipy.sparse as sprs
import scipy.spatial as sptl
from openpnm.core import Base, Workspace, ModelsMixin, logging
from openpnm import topotools
logger = logging.getLogger()
ws = Workspace()


class GenericNetwork(Base, ModelsMixin):
    r"""
    GenericNetwork - Base class to construct pore networks

    """
    def __init__(self, simulation=None, settings={}, **kwargs):
        self.settings.setdefault('prefix', 'net')
        self.settings.update(settings)
        super().__init__(simulation=simulation, **kwargs)
        self._gen_ids()

        # Initialize adjacency and incidence matrix dictionaries
        self._im = {}
        self._am = {}

    def __setitem__(self, prop, value):
        if prop == 'throat.conns':
            if sp.shape(value)[1] != 2:
                logger.error('Wrong size for throat conns!')
            else:
                if sp.any(value[:, 0] > value[:, 1]):
                    logger.warning('Converting throat.conns to be upper \
                                    triangular')
                    value = sp.sort(value, axis=1)
        super().__setitem__(prop, value)

    def __getitem__(self, key):
        if key.split('.')[-1] == self.name:
            element = key.split('.')[0]
            return self[element+'.all']
        if key not in self.keys():
            logger.debug(key + ' not on Network, check on Geometries')
            geoms = self.simulation.geometries().values()
            return self._interleave_data(key, geoms)
        else:
            return super().__getitem__(key)

    def _gen_ids(self):
        if 'pore._id' not in self.keys():
            self['pore._id'] = [str(uuid.uuid4()) for i in self.Ps]
        else:
            # If ids are missing it will from the end of the array...hopefully
            if self['pore._id'][-1] == '':
                inds = sp.where(self['pore._id'] == '')[0]
                temp = [str(uuid.uuid4()) for i in range(len(inds))]
                self['pore._id'][inds] = temp
        if 'throat._id' not in self.keys():
            self['throat._id'] = [str(uuid.uuid4()) for i in self.Ts]
        else:
            # If ids are missing it will from the end of the array...hopefully
            if self['throat._id'][-1] == '':
                inds = sp.where(self['throat._id'] == '')[0]
                temp = [str(uuid.uuid4()) for i in range(len(inds))]
                self['throat._id'][inds] = temp

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
            logger.info('Desired adjacency matrix already present')
            am = self._am[fmt]
        elif self._am.keys():
            logger.info('Desired adjacency matrix not present, converting...')
            am = self._am[list(self._am.keys())[0]]
            tofmt = getattr(am, 'to'+fmt)
            am = tofmt()
            self._am[fmt] = am
        else:
            logger.info('No Adjacency Matrix not present, building...')
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
            An Nt-long array containing the throat values to enter into the
            matrix (in graph theory these are known as the 'weights').  If
            omitted, ones are used to create a standard adjacency matrix
            representing connectivity only.

        fmt : string, optional
            The sparse storage format to return.  Options are:

            **'coo'** : (default) This is the native format of OpenPNM data

            **'lil'** : Enables row-wise slice of the matrix

            **'csr'** : Favored by most linear algebra routines

            **'dok'** : Enables subscript access of locations

        triu : boolean (default is ``False``)
            If ``True``, the returned sparse matrix only contains the upper-
            triangular elements.

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
        >>> temp = pn.create_adjacency_matrix(weights=weights, fmt='csr')

        """
        # Check if provided data is valid
        if weights is None:
            weights = sp.ones((self.Nt,), dtype=int)
        elif sp.shape(weights)[0] != self.Nt:
            raise Exception('Received weights are incorrect length')

        # Append row & col to each other, and data to itself
        conn = self['throat.conns']
        row = conn[:, 0]
        col = conn[:, 1]
        if not triu:
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
        >>> temp = pn.create_incidence_matrix(weights=weights, sprsfmt='csr')
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
            If flatten is True (default) a 1D array of unique pore numbers
            is returned. If flatten is False each location in the the returned
            array contains a sub-arras of neighboring pores for each input
            throat, in the order they were sent.

        mode : string, optional
            Specifies which neighbors should be returned.  The options are:

            **'union'** : (default) All neighbors of the input pores

            **'intersection'** : Only neighbors shared by all input pores

            **'exclusive_or'** : Only neighbors not shared by any input
            pores

        Returns
        -------
        1D array (if flatten is True) or ndarray of arrays (if flatten is
        False)

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
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
        am = self.get_adjacency_matrix(fmt='coo')
        pores = topotools.find_connected_sites(bonds=Ts, am=am,
                                               flatten=flatten, logic=mode)
        return pores

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
        Returns a list the same length as P1 (and P2) with the each element
        containing the throat number that connects the corresponding pores,
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
                            excl_self=True):
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

            **'union'** : (default) All neighbors of the input pores

            **'intersection'** : Only neighbors shared by all input pores

            **'exclusive_or'** : Only neighbors not shared by any input
            pores

        Returns
        -------
        If ``flatten`` is ``True``, returns a 1D array of pore indices filtered
        according to the specified mode.  If ``flatten`` is ``False``, returns
        a list of lists, where each list contains the neighbors of the
        corresponding input pores.

        Notes
        -----
        If ``flatten`` is ``False``, then ``mode`` and ``excl_self`` are
        ignored.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.find_neighbor_pores(pores=[0, 2])
        array([ 1,  3,  5,  7, 25, 27])
        >>> pn.find_neighbor_pores(pores=[0, 1])
        array([ 2,  5,  6, 25, 26])
        >>> pn.find_neighbor_pores(pores=[0, 1], mode='union', excl_self=False)
        array([ 0,  1,  2,  5,  6, 25, 26])
        >>> pn.find_neighbor_pores(pores=[0, 2], flatten=False)
        [[ 1,  5, 25], [ 1,  3,  7, 27]]
        >>> pn.find_neighbor_pores(pores=[0, 2], mode='intersection')
        array([1])
        >>> pn.find_neighbor_pores(pores=[0, 2], mode='exclusive_or')
        array([ 3,  5,  7, 25, 27])
        """
        pores = self._parse_indices(pores)
        neighbors = topotools.find_neighbor_sites(sites=pores, logic=mode,
                                                  am=self.am, flatten=flatten,
                                                  exclude_input=excl_self)
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

            **'union'** : (default) All neighbors of the input pores

            **'intersection'** : Only neighbors shared by all input pores

            **'exclusive_or'** : Only neighbors not shared by any input
            pores

        Returns
        -------
        If ``flatten`` is ``True``, returns a 1D array of throat indices
        filtered according to the specified mode.  If ``flatten`` is ``False``,
        returns a list of lists, where each list contains the neighbors of the
        corresponding input pores.

        Notes
        -----
        If ``flatten`` is ``False``, then ``mode`` and ``excl_self`` are
        ignored.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.find_neighbor_throats(pores=[0, 1])
        array([0, 1, 2, 3, 4, 5])
        >>> pn.find_neighbor_throats(pores=[0, 1], flatten=False)
        array([array([0, 1, 2]), array([0, 3, 4, 5])], dtype=object)
        """
        pores = self._parse_indices(pores)
        neighbors = topotools.find_neighbor_bonds(sites=pores, logic=mode,
                                                  im=self.im, flatten=flatten)
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

            **'exclusive_or'** : The number of neighboring pores that are
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
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.num_neighbors(pores=[0, 1], flatten=False)
        array([3, 4])
        >>> pn.num_neighbors(pores=[0, 2], flatten=True)
        6
        >>> pn.num_neighbors(pores=[0, 2], mode='intersection', flatten=True)
        1
        """
        pores = self._parse_indices(pores)
        # Count number of neighbors
        num = self._find_neighbors(pores, element=element, flatten=flatten,
                                   mode=mode)
        if flatten:
            num = sp.size(num)
        else:
            num = sp.array([sp.size(i) for i in num], dtype=int)
        return num

    def find_nearby_pores(self, pores, r, flatten=False, excl_self=True):
        r"""
        Find all pores within a given radial distance of the input pore(s)
        regardless of whether or not they are toplogically connected.

        Parameters
        ----------
        pores : array_like
            The list of pores for whom nearby neighbors are to be found

        r : scalar
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
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[3, 3, 3])
        >>> pn.find_nearby_pores(pores=[0, 1], r=1)
        array([array([ 1,  5, 25]), array([ 0,  2,  6, 26])], dtype=object)
        >>> pn.find_nearby_pores(pores=[0, 1], r=0.5)
        array([], shape=(2, 0), dtype=int64)
        >>> pn.find_nearby_pores(pores=[0, 1], r=1, flatten=True)
        array([ 2,  5,  6, 25, 26])
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
        [Pn[i].sort() for i in range(0, sp.size(pores))]
        if flatten:  # Convert list of lists to a flat nd-array
            temp = sp.hstack(Pn)
            Pn = sp.unique(Pn)
            if excl_self:  # Remove inputs if necessary
                Pn = Pn[~sp.in1d(Pn, pores)]
        else:  # Convert list of lists to an nd-array of nd-arrays
            if excl_self:  # Remove inputs if necessary
                [Pn[i].remove(pores[i]) for i in range(0, sp.size(pores))]
            temp = [sp.array(Pn[i]) for i in range(0, sp.size(pores))]
            Pn = sp.array(temp)
        if Pn.dtype == float:
            Pn = Pn.astype(sp.int64)
        return Pn
