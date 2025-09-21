import logging
import numpy as np
import scipy.sparse as sprs
import scipy.spatial as sptl
from openpnm.core import Domain
from openpnm import topotools
from openpnm.utils import Docorator
from openpnm.utils import Workspace
import openpnm.models.network as mods
logger = logging.getLogger(__name__)
ws = Workspace()
docstr = Docorator()


@docstr.get_sections(base='NetworkSettings', sections=['Parameters'])
@docstr.dedent
class NetworkSettings:
    r"""

    Parameters
    ----------
    %(BaseSettings.parameters)s
    """


@docstr.get_sections(base='Network', sections=['Parameters'])
@docstr.dedent
class Network(Domain):
    r"""
    This class contains the main functionality used by all networks.

    Parameters
    ----------
    settings : dataclass-like or dict, optional
        User defined settings for the object to override defaults. Can be a
        dataclass-type object with settings stored as attributes or a python
        dicionary of key-value pairs. Settings are stored in the ``settings``
        attribute of the object.
    name : str, optional
        A unique name to assign to the object for easier identification. If
        not given one will be generated.

    coords : array_like (optional)
        An Np-by-3 array of [x, y, z] coordinates for each pore
    conns : array_like (optional)
        An Nt-by-2 array of [head, tail] connections between pores

    Examples
    --------
    >>> import openpnm as op
    >>> import numpy as np

    Create some pore coordinates and connections manually and assign to a
    Network instance. Consider a linear network of 4 pores and 3
    throats:

    ::

        0 ―― 1 ―― 3 ―― 2

    >>> coords = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]])
    >>> conns = np.array([[0, 1], [1, 3], [2, 3]])
    >>> pn = op.network.Network(conns=conns, coords=coords)

    Networks have two required properties: 'pore.coords' and
    'throat.conns'. These arrays indicate the spatial location of each
    pore, and which pores are connected to which. Without these the
    network object cannot function.

    >>> print(pn.props())
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    1     : pore.coords
    2     : throat.conns
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――

    The ``Network`` class has several methods for querying the topology.

    >>> Ps = pn.find_neighbor_pores(pores=1)
    >>> print(Ps)
    [0 3]
    >>> Ts = pn.find_neighbor_throats(pores=[0, 1])
    >>> print(Ts)
    [0 1]
    >>> print(pn.num_neighbors(2))
    [1]

    All of the topological queries are accomplished by inspecting the
    adjacency and incidence matrices. They are created on demand, and are
    stored for future use to save construction time.

    """

    def __init__(self, conns=None, coords=None, name='net', **kwargs):
        super().__init__(name=name, **kwargs)
        self.settings._update(NetworkSettings())
        self._am = {}
        self._im = {}

        if coords is not None:
            coords = np.array(coords)
            self['pore.coords'] = coords

        if conns is not None:
            conns = np.array(conns)
            self['throat.conns'] = conns

        self.add_model(propname='pore.coordination_number',
                       model=mods.coordination_number,
                       regen_mode='deferred',
                       domain='all')
        self.add_model(propname='throat.spacing',
                       model=mods.pore_to_pore_distance,
                       regen_mode='deferred',
                       domain='all')

    def _post_init(self, *args, **kwargs):
        # A collection of functions to run after initialization of a network
        topotools.label_faces(self, label=None)
        drop = ['pore.left', 'pore.right', 'pore.front', 'pore.back',
                'pore.bottom', 'pore.top']
        [self.pop(item) for item in list(self.keys()) if item in drop]

    def __setitem__(self, key, value):
        if key == 'throat.conns':
            if np.any(value[:, 0] > value[:, 1]):
                logger.warning('Converting throat.conns to be upper triangular')
                value = np.sort(value, axis=1)
        super().__setitem__(key, value)

    def get_adjacency_matrix(self, fmt='coo'):
        r"""
        Adjacency matrix in the specified sparse format, with throat IDs
        indicating the non-zero values.

        Parameters
        ----------
        fmt : str, default is 'coo'
            The sparse storage format to return. Options are:
                **'coo'** : This is the native format of OpenPNM's data
                **'lil'** : Enables row-wise slice of the matrix
                **'csr'** : Favored by most linear algebra routines
                **'dok'** : Enables subscript access of locations

        Notes
        -----
        This method will only create the requested matrix in the specified
        format if one is not already saved on the object.  If not present,
        this method will create and return the matrix, as well as store it
        for future use.

        To obtain a matrix with weights other than throat IDs at each
        non-zero location use ``create_adjacency_matrix``.

        To obtain the non-directed graph, with only upper-triangular
        entries, use ``sp.sparse.triu(am, k=1)``.

        """
        # Retrieve existing matrix if available
        if fmt in self._am.keys():
            am = self._am[fmt]
        else:
            am = self.create_adjacency_matrix(weights=self.Ts, fmt=fmt)
            self._am[fmt] = am
        return am

    def get_incidence_matrix(self, fmt='coo'):
        r"""
        Incidence matrix in the specified sparse format, with pore IDs
        indicating the non-zero values.

        Parameters
        ----------
        fmt : str, default is 'coo'
            The sparse storage format to return. Options are:
                **'coo'** : This is the native format of OpenPNM's data
                **'lil'** : Enables row-wise slice of the matrix
                **'csr'** : Favored by most linear algebra routines
                **'dok'** : Enables subscript access of locations

        Notes
        -----
        This method will only create the requested matrix in the specified
        format if one is not already saved on the object. If not present,
        this method will create and return the matrix, as well as store it
        for future use.

        To obtain a matrix with weights other than pore IDs at each
        non-zero location use ``create_incidence_matrix``.

        """
        if fmt in self._im.keys():
            im = self._im[fmt]
        elif self._im.keys():
            im = self._im[list(self._im.keys())[0]]
            tofmt = getattr(im, f"to{fmt}")
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

            If the array is Nt-long, it implies that the matrix is
            symmetric, so the upper and lower triangular regions are
            mirror images. If it is 2*Nt-long then it is assumed that the
            first Nt elements are for the upper triangle, and the last Nt
            element are for the lower triangular.

            If omitted, ones are used to create a standard adjacency
            matrix representing connectivity only.
        fmt : str, optional
            The sparse storage format to return. Options are:
                **'coo'** : This is the native format of OpenPNM's data
                **'lil'** : Enables row-wise slice of the matrix
                **'csr'** : Favored by most linear algebra routines
                **'dok'** : Enables subscript access of locations
        triu : bool, default is ``False``
            If ``True``, the returned sparse matrix only contains the
            upper-triangular elements. This argument is ignored if the
            ``weights`` array is 2*Nt-long.
        drop_zeros : bool, default is ``False``
            If ``True``, applies the ``eliminate_zeros`` method of the
            sparse array to remove all zero locations.

        Returns
        -------
        An adjacency matrix in the specified Scipy sparse format.

        Notes
        -----
        The adjacency matrix is used by OpenPNM for finding the pores
        connected to a give pore or set of pores. Specifically, an
        adjacency matrix has Np rows and Np columns. Each row represents
        a pore, containing non-zero values at the locations corresponding
        to the indices of the pores connected to that pore. The
        ``weights`` argument indicates what value to place at each
        location, with the default being 1's to simply indicate
        connections. Another useful option is throat indices, such that
        the data values on each row indicate which throats are connected
        to the pore.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> weights = np.random.rand(pn.num_throats(), ) < 0.5
        >>> am = pn.create_adjacency_matrix(weights=weights, fmt='csr')

        """
        allowed_weights = [(self.Nt,), (2 * self.Nt,), (self.Nt, 2)]
        # Check if provided data is valid
        if weights is None:
            weights = np.ones((self.Nt,), dtype=int)
        elif np.shape(weights) not in allowed_weights:
            raise Exception('Received weights are of incorrect length')
        weights = np.array(weights)

        # Append row & col to each other, and data to itself
        conn = self['throat.conns']
        row = conn[:, 0]
        col = conn[:, 1]
        if weights.shape == (2 * self.Nt):
            # The flip is necessary since we want [conns.T, reverse(conns).T].T
            row = np.append(row, conn[:, 1])
            col = np.append(col, conn[:, 0])
        elif weights.shape == (self.Nt, 2):
            # The flip is necessary since we want [conns.T, reverse(conns).T].T
            row = np.append(row, conn[:, 1])
            col = np.append(col, conn[:, 0])
            weights = weights.flatten(order='F')
        elif not triu:
            # The flip is necessary since we want [conns.T, reverse(conns).T].T
            row = np.append(row, conn[:, 1])
            col = np.append(col, conn[:, 0])
            weights = np.append(weights, weights)

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
            An array containing the throat values to enter into the matrix
            (in graph theory these are known as the 'weights'). If
            omitted, ones are used to create a standard incidence matrix
            representing connectivity only.
        fmt : str, default is 'coo'
            The sparse storage format to return. Options are:
                **'coo'** : This is the native format of OpenPNM's data
                **'lil'** : Enables row-wise slice of the matrix
                **'csr'** : Favored by most linear algebra routines
                **'dok'** : Enables subscript access of locations
        drop_zeros : bool, default is ``False``
            If ``True``, applies the ``eliminate_zeros`` method of the sparse
            array to remove all zero locations.

        Returns
        -------
        sparse_array
            An incidence matrix in the specified sparse format

        Notes
        -----
        The incidence matrix is a cousin to the adjacency matrix, and used
        by OpenPNM for finding the throats connected to a give pore or set
        of pores. Specifically, an incidence matrix has Np rows and Nt
        columns, and each row represents a pore, containing non-zero
        values at the locations corresponding to the indices of the
        throats connected to that pore. The ``weights`` argument indicates
        what value to place at each location, with the default being 1's
        to simply indicate connections. Another useful option is throat
        indices, such that the data values on each row indicate which
        throats are connected to the pore, though this is redundant as it
        is identical to the locations of non-zeros.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> weights = np.random.rand(pn.num_throats(), ) < 0.5
        >>> im = pn.create_incidence_matrix(weights=weights, fmt='csr')

        """
        # Check if provided data is valid
        if weights is None:
            weights = np.ones((self.Nt,), dtype=int)
        if np.shape(weights)[0] == self.Nt:
            weights = np.append(weights, weights)
        elif np.shape(weights)[0] != 2*self.Nt:
            raise Exception('Received dataset of incorrect length')

        conn = self['throat.conns']
        row = conn[:, 1]
        row = np.append(row, conn[:, 0])
        col = np.arange(self.Nt)
        col = np.append(col, col)

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

    def find_connected_pores(self, throats=[], flatten=False, mode='or'):
        r"""
        Return a list of pores connected to the given list of throats

        Parameters
        ----------
        throats : array_like
            List of throats numbers
        flatten : bool, optional
            If ``True`` (default) a 1D array of unique pore numbers is
            returned. If ``False`` each location in the the returned
            array contains a sub-arras of neighboring pores for each input
            throat, in the order they were sent.
        mode : str
            Specifies logic to filter the resulting list. Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'or'         All neighbors of the input pores. Also accepts 'any'
                         and 'union'.
            'xor'        Only neighbors of one and only one input pore. This
                         is useful for counting the pores that are not shared
                         by any of the input pores. Also accepts
                         'exclusive_or'.
            'xnor'       Neighbors that are shared by two or more input pores.
                         This is equivalent to counting all neighbors
                         with 'or', minus those found with 'xor', and is
                         useful for finding neighbors that the inputs have
                         in common.
            'and'        Only neighbors shared by all input pores. Also
                         accepts 'intersection' and 'all'
            ===========  =====================================================

        Returns
        -------
        1D array (if ``flatten`` is ``True``) or ndarray of arrays (if
        ``flatten`` is ``False``)

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> Ps = pn.find_connected_pores(throats=[0, 1])
        >>> print(Ps)
        [[0 1]
         [1 2]]
        >>> Ps = pn.find_connected_pores(throats=[0, 1], flatten=True)
        >>> print(Ps)
        [0 1 2]

        """
        Ts = self._parse_indices(throats)
        pores = topotools.find_connected_sites(bonds=Ts, network=self,
                                               flatten=flatten, logic=mode)
        return pores

    def find_connecting_throat(self, P1, P2):
        r"""
        Return the throat index connecting pairs of pores.

        Parameters
        ----------
        P1 : array_like
            The indices of the pores whose throats are sought. These can
            be vectors of indices, but must be the same length
        P2 : array_like
            The indices of the pores whose throats are sought. These can
            be vectors of indices, but must be the same length

        Returns
        -------
        list
            Returns a list the same length as P1 (and P2) with the each
            element containing the throat index that connects the
            corresponding pores, or `None`` if pores are not connected.

        Notes
        -----
        The returned list can be converted to an ndarray, which will
        convert the ``None`` values to ``nan``. These can then be found
        using ``numpy.isnan``.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> Ts = pn.find_connecting_throat([0, 1, 2], [2, 2, 2])
        >>> print(Ts)
        [nan  1. nan]

        """
        sites = np.vstack((P1, P2)).T
        Ts = topotools.find_connecting_bonds(sites=sites, network=self)
        return Ts

    def find_neighbor_pores(self, pores, mode='or', flatten=True,
                            include_input=False, asmask=False):
        r"""
        Returns a list of pores that are direct neighbors to the given pore(s)

        Parameters
        ----------
        pores : array_like
            Indices of the pores whose neighbors are sought
        flatten : bool
            If ``True`` (default) the returned result is a compressed
            array of all neighbors.  If ``False``, a list of lists with
            each sub-list containing the neighbors for each input site.
            Note that an *unflattened* list might be slow to generate
            since it is a Python ``list`` rather than a Numpy ``array``.
        include_input : bool
            If ``False`` (default) then the input pores are not included
            in the returned list(s). Note that since pores are not
            neighbors of themselves, the neighbors of pore N will not
            include N, even if this flag is ``True``.
        mode : str
            Specifies logic to filter the resulting list. Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'or'         All neighbors of the input pores. Also accepts 'any'
                         and 'union'.
            'xor'        Only neighbors of one and only one input pore. This
                         is useful for counting the pores that are not shared
                         by any of the input pores. Also accepts
                         'exclusive_or'.
            'xnor'       Neighbors that are shared by two or more input pores.
                         This is equivalent to counting all neighbors
                         with 'or', minus those found with 'xor', and is
                         useful for finding neighbors that the inputs have
                         in common.
            'and'        Only neighbors shared by all input pores. Also
                         accepts 'intersection' and 'all'
            ===========  =====================================================

        asmask : boolean
            If ``False`` (default), the returned result is a list of the
            neighboring pores as indices. If ``True``, the returned result is a
            boolean mask. (Useful for labelling)

        Returns
        -------
        If ``flatten`` is ``True``, returns a 1D array of pore indices
        filtered according to the specified mode.  If ``flatten`` is
        ``False``, returns a list of lists, where each list contains the
        neighbors of the corresponding input pores.

        Notes
        -----
        The ``logic`` options are applied to neighboring pores only, thus
        it is not possible to obtain pores that are part of the global set
        but not neighbors. This is because (a) the list of global pores
        might be very large, and (b) it is not possible to return a list
        of neighbors for each input pores if global pores are considered.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> Ps = pn.find_neighbor_pores(pores=[0, 2])
        >>> print(Ps)
        [ 1  3  5  7 25 27]
        >>> Ps = pn.find_neighbor_pores(pores=[0, 1])
        >>> print(Ps)
        [ 2  5  6 25 26]
        >>> Ps = pn.find_neighbor_pores(pores=[0, 1], mode='or',
        ...                             include_input=True)
        >>> print(Ps)
        [ 0  1  2  5  6 25 26]
        >>> Ps = pn.find_neighbor_pores(pores=[0, 2], flatten=False)
        >>> print(Ps[0])
        [ 1  5 25]
        >>> print(Ps[1])
        [ 1  3  7 27]
        >>> Ps = pn.find_neighbor_pores(pores=[0, 2], mode='xnor')
        >>> print(Ps)
        [1]
        >>> Ps = pn.find_neighbor_pores(pores=[0, 2], mode='xor')
        >>> print(Ps)
        [ 3  5  7 25 27]

        """
        pores = self._parse_indices(pores)
        if np.size(pores) == 0:
            return np.array([], ndmin=1, dtype=int)
        neighbors = topotools.find_neighbor_sites(sites=pores, logic=mode,
                                                  network=self,
                                                  flatten=flatten,
                                                  include_input=include_input)
        if asmask is False:
            return neighbors
        elif flatten is True:
            neighbors = self._tomask(element='pore', indices=neighbors)
            return neighbors
        else:
            raise Exception('Cannot create mask on an unflattened output')

    def find_neighbor_throats(self, pores, mode='or', flatten=True,
                              asmask=False):
        r"""
        Returns a list of throats neighboring the given pore(s)

        Parameters
        ----------
        pores : array_like
            Indices of pores whose neighbors are sought
        flatten : bool, optional
            If ``True`` (default) a 1D array of unique throat indices is
            returned. If ``False`` the returned array contains arrays of
            neighboring throat indices for each input pore, in the order
            they were sent.
        mode : str
            Specifies logic to filter the resulting list. Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'or'         (default) All neighbors of the input throats. Also
                         accepts 'any' and 'union'.
            'xor'        Only neighbors of one and only one input throats. This
                         is useful for counting the pores that are not shared
                         by any of the input pores. Also accepts
                         'exclusive_or'.
            'xnor'       Neighbors that are shared by two or more input
                         throats. This is equivalent to counting all neighbors
                         with 'or', minus those found with 'xor', and is
                         useful for finding neighbors that the inputs have
                         in common.
            'and'        Only neighbors shared by all input throats. Also
                         accepts 'intersection' and 'all'
            ===========  =====================================================

        asmask : boolean
            If ``False`` (default), the returned result is a list of the
            neighboring throats as indices. If ``True``, the returned result is a
            boolean mask. (Useful for labelling)

        Returns
        -------
        If ``flatten`` is ``True``, returns a 1D array of throat indices
        filtered according to the specified mode.  If ``flatten`` is
        ``False``, returns a list of lists, where each list contains the
        neighbors of the corresponding input pores.

        Notes
        -----
        The ``logic`` options are applied to neighboring bonds only, thus
        it is not possible to obtain bonds that are part of the global set
        but not neighbors. This is because (a) the list of global bonds
        might be very large, and (b) it is not possible to return a list
        of neighbors for each input site if global sites are considered.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> Ts = pn.find_neighbor_throats(pores=[0, 1])
        >>> print(Ts)
        [  0   1 100 101 200 201]
        >>> Ts = pn.find_neighbor_throats(pores=[0, 1], flatten=False)
        >>> print(Ts[0])
        [  0 100 200]
        >>> print(Ts[1])
        [  0   1 101 201]

        """
        pores = self._parse_indices(pores)
        if np.size(pores) == 0:
            return np.array([], ndmin=1, dtype=int)
        if flatten is False:
            neighbors = topotools.find_neighbor_bonds(sites=pores, logic=mode,
                                                      network=self,
                                                      flatten=flatten)
        else:
            neighbors = topotools.find_neighbor_bonds(sites=pores, logic=mode,
                                                      network=self, flatten=True)
        if asmask is False:
            return neighbors
        elif flatten is True:
            neighbors = self._tomask(element='throat', indices=neighbors)
            return neighbors
        else:
            raise Exception('Cannot create mask on an unflattened output')

    def num_neighbors(self, pores, mode='or', flatten=False):
        r"""
        Returns the number of neigbhoring pores for each given input pore

        Parameters
        ----------
        pores : array_like
            Pores whose neighbors are to be counted
        flatten : bool, optional
            If ``False`` (default) the number of pores neighboring each
            input pore as an array the same length as ``pores``. If
            ``True`` the sum total number of is counted.
        mode : str
            The logic to apply to the returned count of pores:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'or'         (default) All neighbors of the input pores. Also
                         accepts 'any' and 'union'.
            'xor'        Only neighbors of one and only one input pore. This
                         is useful for counting the pores that are not shared
                         by any of the input pores. Also accepts
                         'exclusive_or'.
            'xnor'       Neighbors that are shared by two or more input pores.
                         This is equivalent to counting all neighbors
                         with 'or', minus those found with 'xor', and is
                         useful for finding neighbors that the inputs have
                         in common.
            'and'        Only neighbors shared by all input pores. Also
                         accepts 'intersection' and 'all'
            ===========  =====================================================

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
        >>> Np = pn.num_neighbors(pores=[0, 1], flatten=False)
        >>> print(Np)
        [3 4]
        >>> Np = pn.num_neighbors(pores=[0, 2], flatten=True)
        >>> print(Np)
        6
        >>> Np = pn.num_neighbors(pores=[0, 2], mode='and', flatten=True)
        >>> print(Np)
        1

        """
        pores = self._parse_indices(pores)
        if flatten:
            # Count number of neighbors
            num = self.find_neighbor_pores(pores, flatten=flatten,
                                           mode=mode, include_input=True)
            num = np.size(num)
        else:   # Could be done much faster if flatten == False
            am = self.create_adjacency_matrix(fmt="csr")
            num = am[pores].sum(axis=1).A1
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
            Controls whether the input pores should be included in the
            list of pores nearby the *other pores* in the input list.
            So if ``pores=[1, 2]`` and 1 and 2 are within ``r`` of each
            other, then 1 will be included in the returned for pores
            near 2, and vice-versa *if* this argument is ``True``.
            The default is ``False``.
        flatten : bool
            If ``True`` returns a single list of all pores that match the
            criteria, otherwise returns an array containing a sub-array for
            each input pore, where each sub-array contains the pores that
            are nearby to each given input pore.  The default is False.

        Returns
        -------
            A list of pores which are within the given spatial distance.
            If a list of N pores is supplied, then a an N-long list of
            such lists is returned.  The returned lists each contain the
            pore for which the neighbors were sought.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[3, 3, 3])
        >>> Ps = pn.find_nearby_pores(pores=[0, 1], r=1)
        >>> print(Ps[0])
        [3 9]
        >>> print(Ps[1])
        [ 2  4 10]
        >>> Ps = pn.find_nearby_pores(pores=[0, 1], r=0.5)
        >>> print(Ps)
        [array([], dtype=int64), array([], dtype=int64)]
        >>> Ps = pn.find_nearby_pores(pores=[0, 1], r=1, flatten=True)
        >>> print(Ps)
        [ 2  3  4  9 10]

        """
        pores = self._parse_indices(pores)
        # Handle an empty array if given
        if np.size(pores) == 0:
            return np.array([], dtype=np.int64)
        if r <= 0:
            raise Exception('Provided distances should be greater than 0')
        # Create kdTree objects
        kd = sptl.cKDTree(self['pore.coords'])
        kd_pores = sptl.cKDTree(self['pore.coords'][pores])
        # Perform search
        Ps_within_r = kd_pores.query_ball_tree(kd, r=r)
        # Remove self from each list
        for i, P in enumerate(Ps_within_r):
            Ps_within_r[i].remove(pores[i])
        # Convert to flattened list by default
        temp = np.concatenate((Ps_within_r))
        Pn = np.unique(temp).astype(np.int64)
        # Remove inputs if necessary
        if include_input is False:
            Pn = Pn[~np.in1d(Pn, pores)]
        # Convert list of lists to a list of ndarrays
        if flatten is False:
            if len(Pn) == 0:  # Deal with no nearby neighbors
                Pn = [np.array([], dtype=np.int64) for i in pores]
            else:
                mask = np.zeros(shape=np.amax((Pn.max(), pores.max())) + 1, dtype=bool)
                mask[Pn] = True
                temp = []
                for item in Ps_within_r:
                    temp.append(np.array(item, dtype=np.int64)[mask[item]])
                Pn = temp
        return Pn

    @property
    def info(self):
        hr = '―' * 78
        errmsg = '--------- ! --------- ! ---------\n'
        pad = 45
        lines = hr + '\n'
        # ---
        lines += 'Number of pores | throats: '.ljust(pad)
        lines += f"{self.Np} | {self.Nt} \n"
        # ---
        lines += 'Pore Coordination (min | mean | max): '.ljust(pad)
        try:
            x = self['pore.coordination_number']
            lines += f"{x.min():.3E} | {x.mean():.3E} | {x.max():.3E} \n"
        except KeyError:
            lines += errmsg
        # ---
        lines += 'Pore Diameter (min | mean | max): '.ljust(pad)
        try:
            x = self['pore.diameter']
            lines += f"{x.min():.3E} | {x.mean():.3E} | {x.max():.3E} \n"
        except KeyError:
            lines += errmsg
        # ---
        lines += 'Pore Spacing (min | mean | max): '.ljust(pad)
        try:
            x = self['throat.spacing']
            lines += f"{x.min():.3E} | {x.mean():.3E} | {x.max():.3E} \n"
        except KeyError:
            lines += errmsg
        # ---
        lines += 'Throat Diameter (min | mean | max): '.ljust(pad)
        try:
            x = self['throat.diameter']
            lines += f"{x.min():.3E} | {x.mean():.3E} | {x.max():.3E} \n"
        except KeyError:
            lines += errmsg
        # ---
        lines += 'Throat Length (min | mean | max): '.ljust(pad)
        try:
            x = self['throat.length']
            lines += f"{x.min():.3E} | {x.mean():.3E} | {x.max():.3E} \n"
        except KeyError:
            lines += errmsg
        # ---
        lines += 'Throat Volume (min | mean | max): '.ljust(pad)
        try:
            x = self['throat.volume']
            lines += f"{x.min():.3E} | {x.mean():.3E} | {x.max():.3E} \n"
        except KeyError:
            lines += errmsg
        # ---
        lines += hr
        print(lines)

    @property
    def conns(self):
        r"""Returns the connectivity matrix of the network."""
        return self['throat.conns']

    @property
    def coords(self):
        r"""Returns the list of pore coordinates of the network."""
        return self['pore.coords']
