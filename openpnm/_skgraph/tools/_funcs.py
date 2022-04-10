import numpy as np
import scipy.sparse as sprs
from scipy.spatial import cKDTree
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
from scipy.sparse import csgraph


# Once a function has been stripped of all its OpenPNM dependent code it
# can be added to this list of functions to import
__all__ = [
    'change_prefix',
    'isoutside',
    'rotate_coords',
    'shear_coords',
    'dimensionality',
    'generate_points_on_sphere',
    'generate_points_in_disk',
    'generate_points_on_circle',
    'find_surface_nodes',
    'find_surface_nodes_cubic',
    'site_to_site_distance',
    'hull_centroid',
    'tri_to_am',
    'vor_to_am',
    'conns_to_am',
    'am_to_im',
    'im_to_am',
    'istriu',
    'istril',
    'istriangular',
    'issymmetric',
]


# This list contains functions that are in the file below but not yet ready
__notyet__ = [
    'is_fully_connected',
    'get_domain_length',
    'get_domain_area',
    'get_shape',
    'get_spacing',
]


def change_prefix(g, old_prefix, new_prefix):
    for item in g.keys():
        if item.startswith(old_prefix):
            key = item.split('.')[1:]
            g[new_prefix + '.' + key] = g.pop(item)
    return g


def isoutside(coords, shape):
    r"""
    Identifies sites that lie outside the specified shape

    Parameters
    ----------
    coords : array_like
        The coordinates which are to be checked
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

    Returns
    -------
    An Np-long mask of ``True`` values indicating pores that lie outside the
    domain.

    Notes
    -----
    If the domain is 2D, either a circle or a square, then the z-dimension
    of ``shape`` should be set to 0.

    """
    # Label external pores for trimming below
    if len(shape) == 1:  # Spherical
        # Find external points
        r = np.sqrt(np.sum(coords**2, axis=1))
        Ps = r > shape[0]
    elif len(shape) == 2:  # Cylindrical
        # Find external pores outside radius
        r = np.sqrt(np.sum(coords[:, [0, 1]]**2, axis=1))
        Ps = r > shape[0]
        # Find external pores above and below cylinder
        if shape[1] > 0:
            Ps = Ps + (coords[:, 2] > shape[1])
            Ps = Ps + (coords[:, 2] < 0)
        else:
            pass
    elif len(shape) == 3:  # Rectilinear
        shape = np.array(shape, dtype=float)
        try:
            lo_lim = shape[:, 0]
            hi_lim = shape[:, 1]
        except IndexError:
            lo_lim = np.array([0, 0, 0])
            hi_lim = shape
        Ps1 = np.any(coords > hi_lim, axis=1)
        Ps2 = np.any(coords < lo_lim, axis=1)
        Ps = Ps1 + Ps2
    return Ps


def rotate_coords(coords, a=0, b=0, c=0, R=None):
    r"""
    Rotates coordinates a given amount about each axis

    Parameters
    ----------
    coords : ndarray
        The site coordinates to be transformed.  ``coords`` must be in 3D,
        but a 2D network can be represented by putting 0's in the missing
        dimension.
    a, b, c : scalar, optional
        The amount in degrees to rotate about the x, y, and z-axis,
        respectively.
    R : array_like, optional
        Rotation matrix.  Must be a 3-by-3 matrix since coordinates are
        always in 3D.  If this is given then `a`, `b`, and `c` are ignored.

    Returns
    -------
    coords : ndarray
        A copy of the given ``coords`` is made and returned to the rotation
        does not occur *in place*.

    See Also
    --------
    shear_coords

    """
    coords = np.copy(coords)
    if R is None:
        if a:
            R = np.array([[1, 0, 0],
                          [0, np.cos(np.deg2rad(a)), -np.sin(np.deg2rad(a))],
                          [0, np.sin(np.deg2rad(a)), np.cos(np.deg2rad(a))]])
            coords = np.tensordot(coords, R, axes=(1, 1))
        if b:
            R = np.array([[np.cos(np.deg2rad(b)), 0, -np.sin(np.deg2rad(b))],
                          [0, 1, 0],
                          [np.sin(np.deg2rad(b)), 0, np.cos(np.deg2rad(b))]])
            coords = np.tensordot(coords, R, axes=(1, 1))
        if c:
            R = np.array([[np.cos(np.deg2rad(c)), -np.sin(np.deg2rad(c)), 0],
                          [np.sin(np.deg2rad(c)), np.cos(np.deg2rad(c)), 0],
                          [0, 0, 1]])
            coords = np.tensordot(coords, R, axes=(1, 1))
    else:
        coords = np.tensordot(coords, R, axes=(1, 1))
    return coords


def shear_coords(coords, ay=0, az=0, bx=0, bz=0, cx=0, cy=0, S=None):
    r"""
    Shears the coordinates a given amount about along axis

    Parameters
    ----------
    coords : ndarray
        The coordinates to be transformed
    ay : scalar
        The factor by which to shear along the x-axis as a function of y
    az : scalar
        The factor by which to shear along the x-axis as a function of z
    bx : scalar
        The factor by which to shear along the y-axis as a function of x
    bz : scalar
        The factor by which to shear along the y-axis as a function of z
    cx : scalar
        The factor by which to shear along the z-axis  as a function of x
    cy : scalar
        The factor by which to shear along the z-axis as a function of y
    S : array_like
        The shear matrix.  Must be a 3-by-3 matrix since pore coordinates are
        always in 3D.  If this is given then the other individual arguments
        are ignored.

    Returns
    -------
    coords : ndarray
        The sheared coordinates.  A copy of the supplied coordinates is made
        so that the operation is not performed *in place*.

    See Also
    --------
    rotate_coords

    Notes
    -----
    The shear along the i *th*-axis is given as i\* = i + aj.  This means
    the new i coordinate is the old one plus some linear factor *a* in the
    j *th* direction.

    The values of ``a``, ``b``, and ``c`` are essentially the inverse of the
    slope to be formed by the neighboring layers of sheared pores.  A value of
    0 means no shear, and neighboring points are stacked directly on top of
    each other; a value of 1 means they form a 45 degree diagonal, and so on.

    If ``S`` is given, then is should be of the form:

    ::

        S = [[1 , ay, az],
             [bx, 1 , bz],
             [cx, cy, 1 ]]

        where any of the off-diagonal components can be 0 meaning no shear

    """
    coords = np.copy(coords)
    if S is None:
        S = np.array([[1, ay, az],
                      [bx, 1, bz],
                      [cx, cy, 1]])
    coords = (S@coords.T).T
    return coords


def dimensionality(coords):
    r"""
    Checks the dimensionality of the network

    Parameters
    ----------
    coords : ndarray
        The coordinates of the sites

    Returns
    -------
    dims : list
        A  3-by-1 array containing ``True`` for each axis that contains
        multiple values, indicating that the pores are spatially distributed
        in that dimension.

    """
    eps = np.finfo(float).resolution
    dims_unique = [not np.allclose(xk, xk.mean(), atol=0, rtol=eps) for xk in coords.T]
    return np.array(dims_unique)


def generate_points_on_sphere(n=100, r=1):
    r"""
    Generates points on a the surface of a sphere with the given radius

    Parameters
    ----------
    n : int or [int, int]
        If a single ``int`` is provided then this number of points will be
        generated using the Fibonacci method to make them approximately
        equally spaced.  If a list of 2 ``int``s is given, they are interpreted
        as the number of meridians and parallels to divide the sphere with.
    r : scalar
        The radius of the sphere on which the points should lie

    Returns
    -------
    coords : ndarray
        An array of x, y, z coordinates for the sphere which will be centered
        on [0, 0, 0]

    """
    if isinstance(n, int):
        i = np.arange(n)
        phi = np.pi * (3 - np.sqrt(5))  # golden angle in radians
        y = 1 - (i / float(n - 1))*2  # y goes from 1 to -1
        radius = np.sqrt(1 - y*y)  # radius at y
        theta = phi * i  # golden angle increment
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        # Convert to spherical coords
        r_, q, p = cart2sph(x, y, z).T
        # Scale the radius then convert back to cartesian
        coords = sph2cart(r=r*r_, theta=q, phi=p)
    else:
        nlat = n[0]
        nlon = n[1]
        lat = []
        lon = []
        for i in range(0, 360, max(1, int(360/nlon))):
            for j in range(max(1, int(180/nlat)), 180, max(1, int(180/nlat))):
                lon.append(i)
                lat.append(j)
        lat.extend([0, 180])
        lon.extend([0, 0])
        theta = np.deg2rad(lon) - np.pi
        phi = np.deg2rad(lat) - np.pi/2
        coords = sph2cart(phi=phi, theta=theta, r=r)
    return coords


def generate_points_in_disk(n=100, r=1):
    r"""
    Generates equally spaced points inside a disk
    """
    indices = np.arange(0, n, dtype=float) + 0.5
    r = np.sqrt(indices/n)
    theta = np.pi*(1 + 5**0.5)*indices
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    return np.vstack((x, y)).T


def generate_points_on_circle(n=100, r=1):
    r"""
    Generates ``n`` equally spaced points on a circle in cartesean coordinates
    """
    theta = np.linspace(0, 2*np.pi, n, endpoint=False)
    r = np.ones_like(theta)*r
    x, y, z = cyl2cart(r=r, theta=theta, z=0)
    return np.vstack((x, y)).T


def cart2sph(x, y, z):
    r"""

    Notes
    -----
    Surprizingly (and annoyingly) this is not built into numpy, for reasons
    discussed `here <https://github.com/numpy/numpy/issues/5228>`_.
    """
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    phi = np.arctan2(z, hxy)
    theta = np.arctan2(y, x)
    return np.vstack((r, theta, phi)).T


def sph2cart(r, theta, phi):
    r"""

    Notes
    -----
    Surprizingly (and annoyingly) this is not built into numpy, for reasons
    discussed `here <https://github.com/numpy/numpy/issues/5228>`_.
    """
    rcos_theta = r * np.cos(phi)
    x = rcos_theta * np.cos(theta)
    y = rcos_theta * np.sin(theta)
    z = r * np.sin(phi)
    return np.vstack((x, y, z)).T


def find_surface_nodes_cubic(coords):
    hits = np.zeros(coords.shape[0], dtype=bool)
    dims = dimensionality(coords)
    for d in range(3):
        if dims[d]:
            hi = np.where(coords[:, d] == coords[:, d].max())[0]
            lo = np.where(coords[:, d] == coords[:, d].min())[0]
            hits[hi] = True
            hits[lo] = True
    return hits


def find_surface_nodes(coords):
    r"""
    Identifies nodes on the outer surface of the domain

    Parameters
    ----------
    coords : ndarray
        The coordinates of nodes in the network

    Returns
    -------
    mask : ndarray
        A boolean array of ``True`` values indicating which nodes were found
        on the surfaces.

    Notes
    -----
    This function generates points around the domain the performs a Delaunay
    tesselation between these points and the network nodes.  Any network
    nodes which are connected to a generated point is considered a surface
    node.
    """
    coords = np.copy(coords)
    shift = np.mean(coords, axis=0)
    coords = coords - shift
    tmp = cart2sph(*coords.T)
    hits = np.zeros(coords.shape[0], dtype=bool)
    r = 2*tmp[0].max()
    dims = dimensionality(coords)
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


def site_to_site_distance(coords, sites1=None, sites2=None):
    r"""
    Find the distance between all pores on set 1 to each pore in set 2

    Parameters
    ----------
    coords : ndarray
        The coordinate of the network sites
    sites1 : array_like
        The site indices of the first set
    sites2 : array_Like
        The site indices of the second set.  It's OK if these indices are
        partially or completely duplicating ``site1``.

    Returns
    -------
    dist : array_like
        A distance matrix with ``len(site1)`` rows and ``len(sites2)`` columns.
        The distance between site *i* in ``site1`` and *j* in ``sites2`` is
        located at *(i, j)* and *(j, i)* in the distance matrix.

    Notes
    -----
    This function computes and returns a distance matrix, which is a dense
    matrix of size Np_1 by Np_2, so can get large.  For distances between
    larger sets a KD-tree approach would be better, which is available in
    ``scipy.spatial``.

    """
    from scipy.spatial.distance import cdist
    p1 = np.array(sites1, ndmin=1)
    p2 = np.array(sites2, ndmin=1)
    return cdist(coords[p1], coords[p2])


def hull_centroid(coords):
    r"""
    Computes centroid of the convex hull enclosing the given coordinates.

    Parameters
    ----------
    coords : Np by 3 ndarray
        Coordinates (xyz)

    Returns
    -------
    centroid : array
        A 3 by 1 Numpy array containing coordinates of the centroid.

    """
    dim = [np.unique(coords[:, i]).size != 1 for i in range(3)]
    hull = ConvexHull(coords[:, dim])
    centroid = coords.mean(axis=0)
    centroid[dim] = hull.points[hull.vertices].mean(axis=0)
    return centroid


def cart2cyl(x, y, z):
    theta = np.arctan2(y, x)
    r = np.hypot(x, y)
    return r, theta, z


def cyl2cart(r, theta, z):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y, z


def to_cylindrical(X, Y, Z):
    r = 2*np.sqrt(X**2 + Y**2)
    theta = 2*np.arctan(Y/X)
    z = Z
    return np.vstack((r, theta, z))


def from_cylindrical(r, theta, z):
    X = r*np.cos(theta)
    Y = r*np.sin(theta)
    Z = z
    return np.vstack((X, Y, Z))


def iscoplanar(coords):
    r"""
    Determines if given pores are coplanar with each other

    Parameters
    ----------
    coords : array_like
        List of pore coords to check for coplanarity.  At least 3 pores are
        required.

    Returns
    -------
    results : bool
        A boolean value of whether given points are coplanar (``True``) or
        not (``False``)

    """
    coords = np.array(coords, ndmin=1)
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


def is_fully_connected(conns, pores_BC=None):
    r"""
    Checks whether network is fully connected, i.e. not clustered.

    Parameters
    ----------
    network : GenericNetwork
        The network whose connectivity to check.
    pores_BC : array_like (optional)
        The pore indices of boundary conditions (inlets/outlets).

    Returns
    -------
    bool
        If ``pores_BC`` is not specified, then returns ``True`` only if
        the entire network is connected to the same cluster. If
        ``pores_BC`` is given, then returns ``True`` only if all clusters
        are connected to the given boundary condition pores.
    """
    am = network.get_adjacency_matrix(fmt='lil').copy()
    temp = csgraph.connected_components(am, directed=False)[1]
    is_connected = np.unique(temp).size == 1
    # Ensure all clusters are part of pores, if given
    if not is_connected and pores_BC is not None:
        am.resize(network.Np + 1, network.Np + 1)
        pores_BC = network._parse_indices(pores_BC)
        am.rows[-1] = pores_BC.tolist()
        am.data[-1] = np.arange(network.Nt, network.Nt + len(pores_BC)).tolist()
        temp = csgraph.connected_components(am, directed=False)[1]
        is_connected = np.unique(temp).size == 1
    return is_connected


def get_spacing(network):
    r"""
    Determine the spacing along each axis of a simple cubic network

    Parameters
    ----------
    network : OpenPNM network
        The network for which spacing is desired

    Returns
    -------
    spacing : ndarray
        The spacing along each axis in the form of ``[Lx, Ly, Lz]``, where
        ``L`` is the physical dimension in the implied units (i.e. meters)

    """
    from openpnm.topotools.generators.tools import get_spacing
    d = {'node.coords': network.coords, 'edge.conns': network.conns}
    spc = get_spacing(d)
    return spc


def get_shape(network):
    r"""
    Determine the shape of each axis of a simple cubic network

    Parameters
    ----------
    network : OpenPNM network
        The network for which shape is desired

    Returns
    -------
    shape : ndarray
        The shape along each axis in the form of ``[Nx, Ny, Nz]`` where
        ``N`` is the number of pores

    """
    from openpnm.topotools.generators.tools import get_shape
    d = {'node.coords': network.coords, 'edge.conns': network.conns}
    shp = get_shape(d)
    return shp


def get_domain_area(network, inlets=None, outlets=None):
    r"""
    Determine the cross sectional area relative to the inlets/outlets.

    Parameters
    ----------
    network : GenericNetwork
        The network object containing the pore coordinates

    inlets : array_like
        The pore indices of the inlets.

    putlets : array_Like
        The pore indices of the outlets.

    Returns
    -------
    area : scalar
        The cross sectional area relative to the inlets/outlets.

    """
    if dimensionality(network).sum() != 3:
        raise Exception('The network is not 3D, specify area manually')
    inlets = network.coords[inlets]
    outlets = network.coords[outlets]
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
    network : GenericNetwork
        The network object containing the pore coordinates

    inlets : array_like
        The pore indices of the inlets.

    outlets : array_Like
        The pore indices of the outlets.

    Returns
    -------
    area : scalar
        The domain length relative to the inlets/outlets.

    """
    inlets = network.coords[inlets]
    outlets = network.coords[outlets]
    if not iscoplanar(inlets):
        print('Detected inlet pores are not coplanar')
    if not iscoplanar(outlets):
        print('Detected inlet pores are not coplanar')
    tree = cKDTree(data=inlets)
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
    # Scan through Delaunay triangulation to retrieve pairs
    indices, indptr = tri.vertex_neighbor_vertices
    for k in range(tri.npoints):
        lil.rows[k] = indptr[indices[k]:indices[k + 1]].tolist()
    # Convert to coo format
    lil.data = lil.rows  # Just a dummy array to make things work properly
    coo = lil.tocoo()
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
    # Create adjacency matrix in lil format for quick matrix construction
    N = vor.vertices.shape[0]
    rc = [[], []]
    for ij in vor.ridge_dict.keys():
        row = vor.ridge_dict[ij].copy()
        # Make sure voronoi cell closes upon itself
        row.append(row[0])
        # Add connections to rc list
        rc[0].extend(row[:-1])
        rc[1].extend(row[1:])
    rc = np.vstack(rc).T
    # Make adj mat upper triangular
    rc = np.sort(rc, axis=1)
    # Remove any pairs with ends at infinity (-1)
    keep = ~np.any(rc == -1, axis=1)
    rc = rc[keep]
    data = np.ones_like(rc[:, 0])
    # Build adj mat in COO format
    M = N = np.amax(rc) + 1
    am = sprs.coo_matrix((data, (rc[:, 0], rc[:, 1])), shape=(M, N))
    # Remove diagonal, and convert to csr remove duplicates
    am = sprs.triu(A=am, k=1, format='csr')
    # The convert back to COO and return
    am = am.tocoo()
    return am


def conns_to_am(conns, shape=None, force_triu=True, drop_diag=True,
                drop_dupes=True, drop_negs=True):
    r"""
    Converts a list of connections into a Scipy sparse adjacency matrix

    Parameters
    ----------
    conns : array_like, N x 2
        The list of site-to-site connections
    shape : list, optional
        The shape of the array.  If none is given then it is inferred from the
        maximum value in ``conns`` array.
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
    if force_triu:  # Sort connections to [low, high]
        conns = np.sort(conns, axis=1)
    if drop_negs:  # Remove connections to -1
        keep = ~np.any(conns < 0, axis=1)
        conns = conns[keep]
    if drop_diag:  # Remove connections of [self, self]
        keep = np.where(conns[:, 0] != conns[:, 1])[0]
        conns = conns[keep]
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


def am_to_im(am):
    r"""
    Convert an adjacency matrix into an incidence matrix
    """
    if am.shape[0] != am.shape[1]:
        raise Exception('Adjacency matrices must be square')
    if am.format != 'coo':
        am = am.tocoo(copy=False)
    conn = np.vstack((am.row, am.col)).T
    row = conn[:, 0]
    data = am.data
    col = np.arange(np.size(am.data))
    if istriangular(am):
        row = np.append(row, conn[:, 1])
        data = np.append(data, data)
        col = np.append(col, col)
    im = sprs.coo.coo_matrix((data, (row, col)), (row.max() + 1, col.max() + 1))
    return im


def im_to_am(im):
    r"""
    Convert an incidence matrix into an adjacency matrix
    """
    raise Exception('This function is not implemented yet')


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
