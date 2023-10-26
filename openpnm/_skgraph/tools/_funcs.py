import numpy as np
import scipy.sparse as sprs
from scipy.spatial import KDTree, distance_matrix
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
from scipy.sparse import csgraph
from openpnm._skgraph.tools import generate_points_on_sphere
from openpnm._skgraph.tools import generate_points_on_circle
from openpnm._skgraph.tools import cart2sph, sph2cart, cart2cyl, cyl2cart


# Once a function has been stripped of all its OpenPNM dependent code it
# can be added to this list of functions to import
__all__ = [
    'get_edge_prefix',
    'get_node_prefix',
    'change_prefix',
    'isoutside',
    'dimensionality',
    'find_surface_nodes',
    'find_surface_nodes_cubic',
    'find_coincident_nodes',
    'internode_distance',
    'tri_to_am',
    'vor_to_am',
    'conns_to_am',
    'dict_to_am',
    'dict_to_im',
    'istriu',
    'istril',
    'isgtriu',
    'istriangular',
    'issymmetric',
    'get_cubic_shape',
    'get_cubic_spacing',
    'is_fully_connected',
    'get_domain_length',
    'get_domain_area',
]


def get_edge_prefix(network):
    r"""
    Determines the prefix used for edge arrays from ``<edge_prefix>.conns``

    Parameters
    ----------
    network : dict
        The network dictionary

    Returns
    -------
    edge_prefix : str
        The value of ``<edge_prefix>`` used in ``g``.  This is found by
        scanning ``g.keys()`` until an array ending in ``'.conns'`` is found,
        then returning the prefix.

    Notes
    -----
    This process is surprizingly fast, on the order of nano seconds, so this
    overhead is worth it for the flexibility it provides in array naming.
    However, since all ``dict`` are now sorted in Python, it may be helpful
    to ensure the ``'conns'`` array is near the beginning of the list.
    """
    for item in network.keys():
        if item.endswith('.conns'):
            return item.split('.')[0]


def get_node_prefix(network):
    r"""
    Determines the prefix used for node arrays from ``<edge_prefix>.coords``

    Parameters
    ----------
    network : dict
        The network dictionary

    Returns
    -------
    node_prefix : str
        The value of ``<node_prefix>`` used in ``g``.  This is found by
        scanning ``g.keys()`` until an array ending in ``'.coords'`` is found,
        then returning the prefix.

    Notes
    -----
    This process is surprizingly fast, on the order of nano seconds, so this
    overhead is worth it for the flexibility it provides in array naming.
    However, since all ``dict`` are now sorted in Python, it may be helpful
    to ensure the ``'conns'`` array is near the beginning of the list.
    """
    for item in network.keys():
        if item.endswith('.coords'):
            return item.split('.')[0]


def change_prefix(network, old_prefix, new_prefix):
    r"""
    Changes the prefix used when generating the graph

    Parameters
    ----------
    network : dict
        The network graph
    old_prefix : str
        The current prefix to change, can either be a node or an edge prefix
    new_prefix : str
        The prefix to use instead

    Returns
    -------
    network : dict
        The graph dictionary will arrays assigned to new keys
    """
    for key in list(network.keys()):
        if key.startswith(old_prefix):
            temp = key.split('.', 1)[1]
            network[new_prefix + '.' + temp] = network.pop(key)
    return network


def isoutside(network, shape, rtol=0.0):
    r"""
    Identifies sites that lie outside the specified shape

    Parameters
    ----------
    network : dict
        The network dictionary. For convenience it is also permissible to just
        supply an N-by-D array of coordinates.
    shape : array_like
        The shape of the domain beyond which points are considered "outside".
        The argument is treated as follows:

        ========== ============================================================
        shape      Interpretation
        ========== ============================================================
        [x, y, z]  A 3D cubic domain of dimension x, y and z with the origin at
                   [0, 0, 0].
        [x, y, 0]  A 2D square domain of size x by y with the origin at
                   [0, 0]
        [r, z]     A 3D cylindrical domain of radius r and height z whose
                   central axis starts at [0, 0, 0]
        [r, 0]     A 2D circular domain of radius r centered on [0, 0] and
                   extending upwards
        [r]        A 3D spherical domain of radius r centered on [0, 0, 0]
        ========== ============================================================

    rtol : scalar or array_like, optional
        Controls how far a node must be from the domain boundary to be
        considered outside. It is applied as a fraction of the domain size as
        ``x[i] > (shape[0] + shape[0]*threshold)`` or
        ``y[i] < (0 - shape[1]*threshold)``.  Discrete threshold values
        can be given for each axis by supplying a list the same size as
        ``shape``.

    Returns
    -------
    mask : boolean ndarray
        A boolean array with ``True`` values indicating nodes that lie outside
        the domain.

    Notes
    -----
    If the domain is 2D, either a circle or a square, then the z-dimension
    of ``shape`` should be set to 0.

    """
    try:
        node_prefix = get_node_prefix(network)
        coords = network[node_prefix+'.coords']
    except AttributeError:
        coords = network
    shape = np.array(shape, dtype=float)
    if np.isscalar(rtol):
        tolerance = np.array([rtol]*len(shape))
    else:
        tolerance = np.array(rtol)
    # Label external pores for trimming below
    if len(shape) == 1:  # Spherical
        # Find external points
        R, Q, P = cart2sph(*coords.T)
        thresh = tolerance[0]*shape[0]
        Ps = R > (shape[0] + thresh)
    elif len(shape) == 2:  # Cylindrical
        # Find external pores outside radius
        R, Q, Z = cart2cyl(*coords.T)
        thresh = tolerance[0]*shape[0]
        Ps = R > shape[0]*(1 + thresh)
        # Find external pores above and below cylinder
        if shape[1] > 0:
            thresh = tolerance[1]*shape[1]
            Ps = Ps + (coords[:, 2] > (shape[1] + thresh))
            Ps = Ps + (coords[:, 2] < (0 - thresh))
        else:
            pass
    elif len(shape) == 3:  # Rectilinear
        thresh = tolerance*shape
        Ps1 = np.any(coords > (shape + thresh), axis=1)
        Ps2 = np.any(coords < (0 - thresh), axis=1)
        Ps = Ps1 + Ps2
    return Ps


def dimensionality(network, cache=True):
    r"""
    Checks the dimensionality of the network

    Parameters
    ----------
    network : dict
        The network dictionary
    cache : boolean, optional (default is True)
        If ``False`` then the dimensionality is recalculated even if it has
        already been calculated and stored in the graph dictionary.

    Returns
    -------
    dims : list
        A  3-by-1 array containing ``True`` for each axis that contains
        multiple values, indicating that the pores are spatially distributed
        in that dimension.

    """
    if cache:
        try:
            return network.params["dimensionality"]
        # union of KeyErroa and AttributeError
        except (KeyError, AttributeError):
            pass
    n = get_node_prefix(network)
    coords = network[n+'.coords']
    eps = np.finfo(float).resolution
    dims_unique = \
        [not np.allclose(xk, xk.mean(), atol=0, rtol=eps) for xk in coords.T]
    if cache:
        network["params.dimensionality"] = np.array(dims_unique)
    return np.array(dims_unique)


def find_surface_nodes_cubic(network):
    r"""
    Identifies nodes on the outer surface of the domain assuming a cubic domain
    to save time

    Parameters
    ----------
    network : dict
        The graph dictionary

    Returns
    -------
    mask : ndarray
        A boolean array of ``True`` values indicating which nodes were found
        on the surfaces.
    """
    node_prefix = get_node_prefix(network)
    coords = network[node_prefix+'.coords']
    hits = np.zeros(coords.shape[0], dtype=bool)
    dims = dimensionality(network)
    for d in range(3):
        if dims[d]:
            hi = np.where(coords[:, d] == coords[:, d].max())[0]
            lo = np.where(coords[:, d] == coords[:, d].min())[0]
            hits[hi] = True
            hits[lo] = True
    return hits


def find_coincident_nodes(network):
    r"""
    Finds nodes with identical coordinates

    Parameters
    ----------
    network : dict
        The network dictionary

    Returns
    -------
    duplicates : list of ndarrays
        A list with each sublist indicating the indices of nodes that share
        a common set of coordinates

    Notes
    -----
    This function works by computing a ``hash`` of the coordinates then finding
    all nodes with equivalent hash values. Hashes are supposed to be unique
    but they occassionally "collide", meaning nodes may be identified as
    coincident that are not.
    """
    node_prefix = get_node_prefix(network)
    coords = network[node_prefix+'.coords']
    hashed = [hash(row.tobytes()) for row in coords]
    uniq, counts = np.unique(hashed, return_counts=True)
    hits = np.where(counts > 1)[0]
    dupes = []
    for item in hits:
        dupes.append(np.where(hashed == uniq[item])[0])
    return dupes


def find_surface_nodes(network):
    r"""
    Identifies nodes on the outer surface of the domain using a Delaunay
    tessellation

    Parameters
    ----------
    network : dict
        The network dictionary

    Returns
    -------
    mask : ndarray
        A boolean array of ``True`` values indicating which nodes were found
        on the surfaces.

    Notes
    -----
    This function generates points around the domain the performs a Delaunay
    tesselation between these points and the network nodes.  Any network
    nodes which are connected to a generated points is considered a surface
    node.

    """
    node_prefix = get_node_prefix(network)
    coords = np.copy(network[node_prefix+'.coords'])
    shift = np.mean(coords, axis=0)
    coords = coords - shift
    tmp = cart2sph(*coords.T)
    hits = np.zeros(coords.shape[0], dtype=bool)
    r = 2*tmp[0].max()
    dims = dimensionality(network)
    if sum(dims) == 1:
        hi = np.where(coords[:, dims] == coords[:, dims].max())[0]
        lo = np.where(coords[:, dims] == coords[:, dims].min())[0]
        hits[hi] = True
        hits[lo] = True
        return hits
    if sum(dims) == 2:
        markers = generate_points_on_circle(n=max(10, int(coords.shape[0]/10)), r=r)
        pts = np.vstack((coords[:, dims], markers))
    else:
        markers = generate_points_on_sphere(n=max(10, int(coords.shape[0]/10)), r=r)
        pts = np.vstack((coords, markers))
    tri = Delaunay(pts, incremental=False)
    (indices, indptr) = tri.vertex_neighbor_vertices
    for k in range(coords.shape[0], tri.npoints):
        neighbors = indptr[indices[k]:indices[k+1]]
        inds = np.where(neighbors < coords.shape[0])
        hits[neighbors[inds]] = True
    return hits


def internode_distance(network, inds_1=None, inds_2=None):
    r"""
    Find the distance between all nodes on set 1 to each node in set 2

    Parameters
    ----------
    network : dict
        The network dictionary
    inds_1 : array_like
        A list containing the indices of the first set of nodes
    inds_2 : array_Like
        A list containing the indices of the first set of nodes.  It's OK if
        these indices are partially or completely duplicating ``nodes1``.

    Returns
    -------
    dist : array_like
        A distance matrix with ``len(site1)`` rows and ``len(sites2)`` columns.
        The distance between site *i* in ``site1`` and *j* in ``sites2`` is
        located at *(i, j)* and *(j, i)* in the distance matrix.

    Notes
    -----
    This function computes and returns a distance matrix, so can get large.
    For distances between larger sets a KD-tree approach would be better,
    which is available in ``scipy.spatial``.

    """
    node_prefix = get_node_prefix(network)
    coords = network[node_prefix+'.coords']
    p1 = np.array(inds_1, ndmin=1)
    p2 = np.array(inds_2, ndmin=1)
    return distance_matrix(coords[p1], coords[p2])


def iscoplanar(network):
    r"""
    Determines if specified nodes are coplanar with each other

    Parameters
    ----------
    network : dict
        The graph dictionary

    Returns
    -------
    flag : bool
        A boolean value of whether given nodes are coplanar (``True``) or
        not (``False``)

    """
    node_prefix = get_node_prefix(network)
    coords = network[node_prefix + '.coords']
    if np.shape(coords)[0] < 3:
        raise Exception('At least 3 input pores are required')

    Px = coords[:, 0]
    Py = coords[:, 1]
    Pz = coords[:, 2]

    # Do easy check first, for common coordinate
    if np.shape(np.unique(Px))[0] == 1:
        return True
    if np.shape(np.unique(Py))[0] == 1:
        return True
    if np.shape(np.unique(Pz))[0] == 1:
        return True

    # Perform rigorous check using vector algebra
    # Grab first basis vector from list of coords
    n1 = np.array((Px[1] - Px[0], Py[1] - Py[0], Pz[1] - Pz[0])).T
    n = np.array([0.0, 0.0, 0.0])
    i = 1
    while n.sum() == 0:
        if i >= (np.size(Px) - 1):
            return False
        # Chose a secon basis vector
        n2 = np.array((Px[i+1] - Px[i], Py[i+1] - Py[i], Pz[i+1] - Pz[i])).T
        # Find their cross product
        n = np.cross(n1, n2)
        i += 1
    # Create vectors between all other pairs of points
    r = np.array((Px[1:-1] - Px[0], Py[1:-1] - Py[0], Pz[1:-1] - Pz[0]))
    # Ensure they all lie on the same plane
    n_dot = np.dot(n, r)

    return bool(np.sum(np.absolute(n_dot)) == 0)


def is_fully_connected(network, inds=None):
    r"""
    Checks whether graph is fully connected, i.e. not clustered

    Parameters
    ----------
    network : dict
        The network dictionary
    inds : array_like (optional)
        The indices of boundary nodes (i.e. inlets/outlets). If this is given
        the multiple sample spanning clusters will count as fully connected.

    Returns
    -------
    flag : bool
        If ``inds`` is not specified, then returns ``True`` only if
        the entire network is connected to the same cluster. If
        ``inds`` is given, then returns ``True`` only if all clusters
        are connected to the given boundary nodes.
    """
    am = dict_to_am(network)
    am = am.tolil()
    inds = np.array(inds)
    temp = csgraph.connected_components(am, directed=False)[1]
    is_connected = np.unique(temp).size == 1
    Np = am.shape[0]
    Nt = int(am.nnz/2)
    # Ensure all clusters are part of inds, if given
    if not is_connected and inds is not None:
        am.resize(Np + 1, Np + 1)
        am.rows[-1] = inds.tolist()
        am.data[-1] = np.arange(Nt, Nt + len(inds)).tolist()
        temp = csgraph.connected_components(am, directed=False)[1]
        is_connected = np.unique(temp).size == 1
    return is_connected


def get_cubic_spacing(network):
    r"""
    Determine spacing of a cubic network

    Parameters
    ----------
    network : dict
        The network dictionary

    Returns
    -------
    spacing : ndarray
        An array containing the spacing between nodes in each direction

    """
    node_prefix = get_node_prefix(network)
    edge_prefix = get_edge_prefix(network)
    coords = network[node_prefix+'.coords']
    conns = network[edge_prefix+'.conns']
    # Find Network spacing
    C12 = coords[conns]
    mag = np.linalg.norm(np.diff(C12, axis=1), axis=2)
    unit_vec = np.around(np.squeeze(np.diff(C12, axis=1)) / mag, decimals=14)
    spacing = [0, 0, 0]
    dims = dimensionality(network)
    # Ensure vectors point in n-dims unique directions
    c = {tuple(row): 1 for row in unit_vec}
    mag = np.atleast_1d(mag.squeeze()).astype(float)
    if len(c.keys()) > sum(dims):
        raise Exception(
            "Spacing is undefined when throats point in more directions"
            " than network has dimensions."
        )
    for ax in [0, 1, 2]:
        if dims[ax]:
            inds = np.where(unit_vec[:, ax] == unit_vec[:, ax].max())[0]
            temp = np.unique(mag[inds])
            if not np.allclose(temp, temp[0]):
                raise Exception("A unique value of spacing could not be found.")
            spacing[ax] = temp[0]
    return np.array(spacing)


def get_cubic_shape(network):
    r"""
    Determine shape of a cubic network

    Parameters
    ----------
    network : dict
        The network dictionary

    Returns
    -------
    shape : ndarray
        An array containing the shape of the network each direction

    """
    node_prefix = get_node_prefix(network)
    coords = network[node_prefix+'.coords']
    L = np.ptp(coords, axis=0)
    mask = L.astype(bool)
    S = get_cubic_spacing(network)
    shape = np.array([1, 1, 1], int)
    shape[mask] = L[mask] / S[mask] + 1
    return shape


def get_domain_area(network, inlets=None, outlets=None):
    r"""
    Determine the cross sectional area relative to the inlets/outlets.

    Parameters
    ----------
    network : dict
        The network dictionary
    inlets : array_like
        The indices of the inlets
    outlets : array_Like
        The indices of the outlets

    Returns
    -------
    area : scalar
        The cross sectional area relative to the inlets/outlets.

    """
    if dimensionality(network).sum() != 3:
        raise Exception('The network is not 3D, specify area manually')
    node_prefix = get_node_prefix(network)
    coords = network[node_prefix+'.coords']
    inlets = coords[inlets]
    outlets = coords[outlets]
    if not iscoplanar(inlets):
        print('Detected inlet pores are not coplanar')
    if not iscoplanar(outlets):
        print('Detected outlet pores are not coplanar')
    Nin = np.ptp(inlets, axis=0) > 0
    if Nin.all():
        print('Detected inlets are not oriented along a principle axis')
    Nout = np.ptp(outlets, axis=0) > 0
    if Nout.all():
        print('Detected outlets are not oriented along a principle axis')
    hull_in = ConvexHull(points=inlets[:, Nin])
    hull_out = ConvexHull(points=outlets[:, Nout])
    if hull_in.volume != hull_out.volume:
        print('Inlet and outlet faces are different area')
    area = hull_in.volume  # In 2D: volume=area, area=perimeter
    return area


def get_domain_length(network, inlets=None, outlets=None):
    r"""
    Determine the domain length relative to the inlets/outlets.

    Parameters
    ----------
    network : dict
        The network dictionary
    inlets : array_like
        The pore indices of the inlets.
    outlets : array_Like
        The pore indices of the outlets.

    Returns
    -------
    area : scalar
        The domain length relative to the inlets/outlets.

    """
    node_prefix = get_node_prefix(network)
    coords = network[node_prefix+'.coords']
    inlets = coords[inlets]
    outlets = coords[outlets]
    if not iscoplanar(inlets):
        print('Detected inlet pores are not coplanar')
    if not iscoplanar(outlets):
        print('Detected inlet pores are not coplanar')
    tree = KDTree(data=inlets)
    Ls = np.unique(np.float64(tree.query(x=outlets)[0]))
    if not np.allclose(Ls, Ls[0]):
        print('A unique value of length could not be found')
    length = Ls[0]
    return length


def tri_to_am(tri):
    r"""
    Given a Delaunay triangulation object from Scipy's ``spatial`` module,
    converts to a sparse adjacency matrix network representation.

    Parameters
    ----------
    tri : Delaunay Triangulation Object
        This object is produced by ``scipy.spatial.Delaunay``

    Returns
    -------
    A sparse adjacency matrix in COO format.  The network is undirected
    and unweighted, so the adjacency matrix is upper-triangular and all the
    weights are set to 1.

    """
    # Create an empty list-of-list matrix
    lil = sprs.lil_matrix((tri.npoints, tri.npoints))
    # Scan through Delaunay triangulation object to retrieve pairs
    indices, indptr = tri.vertex_neighbor_vertices
    if 1:  # Original way
        for k in range(tri.npoints):
            lil.rows[k] = indptr[indices[k]:indices[k + 1]].tolist()
        # Convert to coo format
        lil.data = lil.rows  # Just a dummy array to make things work properly
        coo = lil.tocoo()
    else:  # Alternative way, about 10x SLOWER
        coo = [[], []]
        for k in range(tri.npoints):
            # lil.rows[k] = indptr[indices[k]:indices[k + 1]].tolist()
            col = indptr[indices[k]:indices[k+1]]
            coo[0].extend(np.ones_like(col)*k)
            coo[1].extend(col)
        coo = sprs.coo_matrix((np.ones_like(coo[0]), (coo[0], coo[1])))
    # Set weights to 1's
    coo.data = np.ones_like(coo.data)
    # Remove diagonal, and convert to csr remove duplicates
    am = sprs.triu(A=coo, k=1, format='csr')
    # The convert back to COO and return
    am = am.tocoo()
    return am


def vor_to_am(vor):
    r"""
    Given a Voronoi tessellation object from Scipy's ``spatial`` module,
    converts to a sparse adjacency matrix network representation in COO format.

    Parameters
    ----------
    vor : Voronoi Tessellation object
        This object is produced by ``scipy.spatial.Voronoi``

    Returns
    -------
    A sparse adjacency matrix in COO format.  The network is undirected
    and unweighted, so the adjacency matrix is upper-triangular and all the
    weights are set to 1.

    """
    if 0:  # Original way, 2X slower
        rc = [[], []]
        for ij in vor.ridge_dict.keys():
            row = vor.ridge_dict[ij].copy()
            # Make sure voronoi cell closes upon itself
            row.append(row[0])
            # Add connections to rc list
            rc[0].extend(row[:-1])
            rc[1].extend(row[1:])
        rc = np.vstack(rc).T
    else:
        v = vor.ridge_vertices.copy()
        # Add row [0] to close the facet on itself, add -1 to break connection to
        # next facet in list as connections with -1 get deleted
        _ = [row.extend([row[0], -1]) for row in v]
        v = np.hstack(v)
        rc = np.vstack((v[:-1], v[1:])).T
    mask = np.any(rc < 0, axis=1)
    rc = rc[~mask]
    rc = np.sort(rc, axis=1)
    am = conns_to_am(rc)
    return am


def dict_to_am(network, weights=None):
    r"""
    Convert a graph dictionary into a ``scipy.sparse`` adjacency matrix in
    COO format

    Parameters
    ----------
    network : dict
        A network dictionary
    weights : ndarray, optional
        The weight values to use for the connections. If not provided
        then 1's are assumed.

    Returns
    -------
    am : sparse matrix
        The sparse adjacency matrix in COO format

    Notes
    -----
    If the edge connections in ``g`` are in upper-triangular form, then the
    graph is assumed to be undirected and the returned adjacency matrix is
    symmetrical (i.e. the triu entries are reflected in tril). If any edges
    are found in the lower triangle, then the returned adjacency matrix is
    unsymmetrical.

    Multigraphs (i.e. duplicate connections between nodes) are not suported,
    but this is not checked for here to avoid overhead since this function is
    called frequently.

    """
    edge_prefix = get_edge_prefix(network)
    node_prefix = get_node_prefix(network)
    conns = np.copy(network[edge_prefix+'.conns'])
    shape = [network[node_prefix+'.coords'].shape[0]]*2
    if weights is None:
        weights = np.ones_like(conns[:, 0], dtype=int)
    if isgtriu(network):  # If graph is triu, then it is assumed to be undirected
        conns = np.vstack((conns, np.fliplr(conns)))  # Reflect to tril
        data = np.ones_like(conns[:, 0], dtype=int)  # Generate fake data
        am = sprs.coo_matrix((data, (conns[:, 0], conns[:, 1])), shape=shape)
        am.data = np.hstack((weights, weights))
    else:
        am = sprs.coo_matrix((weights, (conns[:, 0], conns[:, 1])), shape=shape)
    return am


def dict_to_im(network):
    r"""
    Convert a graph dictionary into a ``scipy.sparse`` incidence matrix in COO
    format

    Parameters
    ----------
    network : dict
        The network dictionary

    Returns
    -------
    im : sparse matrix
        The sparse incidence matrix in COO format

    Notes
    -----
    Rows correspond to nodes and columns correspond to edges. Each column
    has 2 nonzero values indicating which 2 nodes are connected by the
    corresponding edge. Each row contains an arbitrary number of nonzeros
    whose locations indicate which edges are directly connected to the
    corresponding node.
    """
    edge_prefix = get_edge_prefix(network)
    node_prefix = get_node_prefix(network)
    conns = network[edge_prefix+'.conns']
    coords = network[node_prefix+'.coords']
    if isgtriu(network):
        data = np.ones(2*conns.shape[0], dtype=int)
        shape = (coords.shape[0], conns.shape[0])
        temp = np.arange(conns.shape[0])
        cols = np.vstack((temp, temp)).T.flatten()
        rows = conns.flatten()
        im = sprs.coo_matrix((data, (rows, cols)), shape=shape)
    else:
        raise Exception('This function is not implemented for directed graphs')
    return im


def ismultigraph(network):
    r"""
    Checks if graph contains multiple connections between any pair of nodes

    Parameters
    ----------
    network : dict
        The network dictionary

    Returns
    -------
    flag : bool
        Returns ``True`` if any pair of nodes is connected by more than one
        edge.
    """
    edge_prefix = get_edge_prefix(network)
    node_prefix = get_node_prefix(network)
    conns = network[edge_prefix+'.conns']
    coords = network[node_prefix+'.coords']
    data = np.ones_like(conns[:, 0], dtype=int)
    shape = 2*[coords.shape[0]]
    am = sprs.coo_matrix((data, (conns[:, 0], conns[:, 1])), shape=shape)
    am.sum_duplicates()
    return np.any(am.data > 1)


def isgtriu(network):
    r"""
    Determines if graph connections are in upper triangular format

    Parameters
    ----------
    network : dict
        The network dictionary

    Returns
    -------
    flag : bool
        Returns ``True`` if *all* rows in "conns" are ordered as [lo, hi]
    """
    edge_prefix = get_edge_prefix(network)
    conns = network[edge_prefix+'.conns']
    return np.all(conns[:, 0] < conns[:, 1])


def to_triu(network):
    r"""
    Adjusts conns array to force into upper triangular form

    Parameters
    ----------
    network : dict
        The network dictionary

    Returns
    -------
    network : dict
        The graph dictionary with edge connections updated

    Notes
    -----
    This does not check for the creation of duplicate connections
    """
    edge_prefix = get_edge_prefix(network)
    conns = network[edge_prefix+'.conns']
    network[edge_prefix+'.conns'] = np.sort(conns, axis=1)
    return network


def conns_to_am(conns, shape=None, force_triu=True, drop_diag=True,
                drop_dupes=True, drop_negs=True):
    r"""
    Converts a list of connections into a Scipy sparse adjacency matrix

    Parameters
    ----------
    conns : array_like, N x 2
        The list of site-to-site connections
    shape : list, optional
        The shape of the array.  If none is given then it is taken as 1 + the
        maximum value in ``conns``.
    force_triu : bool
        If True (default), then all connections are assumed undirected, and
        moved to the upper triangular portion of the array
    drop_diag : bool
        If True (default), then connections from a site and itself are removed.
    drop_dupes : bool
        If True (default), then all pairs of sites sharing multiple connections
        are reduced to a single connection.
    drop_negs : bool
        If True (default), then all connections with one or both ends pointing
        to a negative number are removed.

    Returns
    -------
    am : ndarray
        A sparse adjacency matrix in COO format

    """
    if drop_negs:  # Remove connections to -1
        keep = ~np.any(conns < 0, axis=1)
        conns = conns[keep]
    if drop_diag:  # Remove connections of [self, self]
        keep = np.where(conns[:, 0] != conns[:, 1])[0]
        conns = conns[keep]
    if force_triu:  # Sort connections to [low, high]
        conns = np.sort(conns, axis=1)
    # Now convert to actual sparse array in COO format
    data = np.ones_like(conns[:, 0], dtype=int)
    if shape is None:
        N = conns.max() + 1
        shape = (N, N)
    am = sprs.coo_matrix((data, (conns[:, 0], conns[:, 1])), shape=shape)
    if drop_dupes:  # Convert to csr and back too coo
        am = am.tocsr()
        am = am.tocoo()
    # Perform one last check on adjacency matrix
    missing = np.where(np.bincount(conns.flatten()) == 0)[0]
    if np.size(missing) or np.any(am.col.max() < (shape[0] - 1)):
        print('Some nodes are not connected to any bonds')
    return am


def istriu(am):
    r"""
    Returns ``True`` if the sparse adjacency matrix is upper triangular
    """
    if am.shape[0] != am.shape[1]:
        print('Matrix is not square, triangularity is irrelevant')
        return False
    if am.format != 'coo':
        am = am.tocoo(copy=False)
    return np.all(am.row <= am.col)


def istril(am):
    r"""
    Returns ``True`` if the sparse adjacency matrix is lower triangular
    """
    if am.shape[0] != am.shape[1]:
        print('Matrix is not square, triangularity is irrelevant')
        return False
    if am.format != 'coo':
        am = am.tocoo(copy=False)
    return np.all(am.row >= am.col)


def istriangular(am):
    r"""
    Returns ``True`` if the sparse adjacency matrix is either upper or lower
    triangular
    """
    if am.format != 'coo':
        am = am.tocoo(copy=False)
    return istril(am) or istriu(am)


def issymmetric(am):
    r"""
    A method to check if a square matrix is symmetric
    Returns ``True`` if the sparse adjacency matrix is symmetric
    """
    if am.shape[0] != am.shape[1]:
        print('Matrix is not square, symmetrical is irrelevant')
        return False
    if am.format != 'coo':
        am = am.tocoo(copy=False)
    if istril(am) or istriu(am):
        return False
    # Compare am with its transpose, element wise
    sym = ((am != am.T).size) == 0
    return sym
