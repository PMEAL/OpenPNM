"""
Tools
-----

"""
import numpy as np


__all__ = [
    'trim',
    'join',
    'parse_points',
    'get_spacing',
    'get_shape',
    'add_all_label',
    'label_surface_nodes',
    'label_faces',
    'crop',
    'template_sphere_shell',
    'template_cylinder_annulus',
]


def trim(network, bonds=None, sites=None):
    r"""
    Remove given pores or throats from a network

    Parameters
    ----------
    network : dictionary
        A dictionary containing 'vert.coords' and 'edge.conns', describing
        the network
    sites : array_like
        The list of vertcies to trim.
    bonds : array_like
        The list of edges to trim

    Returns
    -------
    network : dict
        The ``network`` object with the specified vertices/edges removed

    Notes
    -----
    Only one of ``vert_ids`` or ``edge_ids`` can be given.  To trim both types,
    call the fuction twice. This function renumbers the ``'edge.conns'``
    array when the vertices being pointed to are removed.

    """
    if (bonds is not None) and (sites is not None):
        raise Exception('Cannot trim pores and throats at the same time')
    if bonds is not None:
        edge_ids = np.atleast_1d(bonds)
        keep = np.ones(network['conns'].shape[0], dtype=bool)
        keep[edge_ids] = False
        for item in network.keys():
            if item.startswith('edge'):
                network[item] = network[item][keep]
    elif sites is not None:
        vert_ids = np.atleast_1d(sites)
        if vert_ids.dtype == bool:
            vert_ids = np.where(vert_ids)[0]
        keep = np.ones(network['coords'].shape[0], dtype=bool)
        keep[vert_ids] = False
        for item in network.keys():
            if item.startswith('vert'):
                network[item] = network[item][keep]
        # Remove edges
        edges = np.any(np.isin(network['conns'], vert_ids), axis=1)
        network = trim(network, bonds=edges)
        # Renumber throat conns
        remapping = np.cumsum(keep) - 1
        network['conns'] = remapping[network['conns']]
    return network


def join(net1, net2, L_max=0.99):
    r"""
    Joins two networks together topologically including new connections

    Parameters
    ----------
    net1 : dictionary
        A dictionary containing 'vert.coords' and 'edge.conns'.
    net2 : dictionary
        A dictionary containing 'vert.coords' and 'edge.conns'
    L_max : float
        The maximum distance between vertices below which they are called
        neighbors

    Returns
    -------
    network : dict
        A dictionary containing 'vert.coords' vertices from both ``net1`` and
        ``net2``, and ``edge.conns`` with original connections plus new ones
        found during the join process.

    Notes
    -----
    This function uses ``scipy.spatial.KDTree``.

    """
    # Perform neighbor query
    from scipy.spatial import KDTree
    t1 = KDTree(net1['vert.coords'])
    t2 = KDTree(net2['vert.coords'])
    pairs = t1.query_ball_tree(t2, r=0.99)
    # Combine existing network data
    net3 = {}
    Np1 = net1['vert.coords'].shape[0]
    Np2 = net2['vert.coords'].shape[0]
    net3['vert.coords'] = np.vstack((net1.pop('vert.coords'),
                                     net2.pop('vert.coords')))
    net3['edge.conns'] = np.vstack((net1.pop('edge.conns'),
                                    net2.pop('edge.conns') + Np1))
    # Convert kdtree result into new connections
    nnz = sum([len(row) for row in pairs])
    conns = np.zeros((nnz, 2), dtype=int)
    i = 0
    for j, row in enumerate(pairs):
        for col in row:
            conns[i, :] = j, col + Np1
            i += 1
    # Add new connections to network
    net3['edge.conns'] = np.vstack((net3.pop('edge.conns'), conns))
    # Finally, expand any other data arrays on given networks
    keys = set(net1.keys()).union(net2.keys())
    for item in keys:
        temp1 = net1.pop(item, np.zeros(Np1)*np.nan)
        temp2 = net2.pop(item, np.zeros(Np2)*np.nan)
        net3[item] = np.concatenate((temp1, temp2), axis=0)
    return net3


def parse_points(shape, points, reflect=False):
    from openpnm.topotools import generate_base_points
    # Deal with input arguments
    if isinstance(points, int):
        points = generate_base_points(num_points=points,
                                      domain_size=shape,
                                      reflect=reflect)
    else:
        # Should we check to ensure that points are reflected?
        points = np.array(points)
    return points


def get_spacing(network):
    r"""
    Determine spacing of a cubic network

    Parameters
    ----------
    network : dictionary
        A network dictionary containing 'vert.coords' and 'edge.conns'

    Returns
    -------
    spacing : ndarray
        An array containing the spacing between vertices in each direction

    Notes
    -----
    This function only works on simple cubic networks with no boundary
    vertices. If a unique spacing cannot be found in each direction,
    and/or the edges are not all oriented perpendicularly, exceptions
    will be raised.

    """
    from openpnm.topotools import dimensionality
    # Find Network spacing
    P12 = network["edge.conns"]
    C12 = network["vert.coords"][P12]
    mag = np.linalg.norm(np.diff(C12, axis=1), axis=2)
    unit_vec = np.around(np.squeeze(np.diff(C12, axis=1)) / mag, decimals=14)
    spacing = [0, 0, 0]
    dims = dimensionality(coords=network['vert.coords'])
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


def get_shape(network):
    L = np.ptp(network["vert.coords"], axis=0)
    mask = L.astype(bool)
    S = get_spacing(network)
    shape = np.array([1, 1, 1], int)
    shape[mask] = L[mask] / S[mask] + 1
    return shape


def add_all_label(network):
    network['pore.all'] = np.ones(network['pore.coords'].shape[0], dtype=bool)
    network['throat.all'] = np.ones(network['throat.conns'].shape[0], dtype=bool)
    return network


def label_surface_nodes(network):
    r"""
    """
    from openpnm.topotools import dimensionality
    hits = np.zeros_like(network.Ps, dtype=bool)
    dims = dimensionality(network)
    mn = np.amin(network["vert.coords"], axis=0)
    mx = np.amax(network["vert.coords"], axis=0)
    for ax in np.where(dims)[0]:
        if dims[ax]:
            hits += network["vert.coords"][:, ax] <= mn[ax]
            hits += network["vert.coords"][:, ax] >= mx[ax]
    network["vert.surface"] = hits
    return network


def label_faces(network, threshold=0.05):
    r"""
    Label the vertices sitting on the faces of the domain in accordance with
    the conventions used for cubic etc.
    """
    from openpnm.topotools import dimensionality
    dims = dimensionality(network)
    coords = np.around(network['vert.coords'], decimals=10)
    min_labels = ['front', 'left', 'bottom']
    max_labels = ['back', 'right', 'top']
    min_coords = np.amin(coords, axis=0)
    max_coords = np.amax(coords, axis=0)
    for ax in np.where(dims)[0]:
        network['vert.' + min_labels[ax]] = coords[:, ax] <= threshold*min_coords[ax]
        network['vert.' + max_labels[ax]] = coords[:, ax] >= (1-threshold)*max_coords[ax]
    return network


def crop(network, shape, mode='full'):
    r"""
    Removes vertices that lie outside the specifed shape

    Parameters
    ----------
    network : dict
        Dictionary containing ``vert.coords`` and ``edge.conns`` arrays
    shape : array_like
        The [x, y, z] shape of the domain beyond which trimming should be
        applied
    mode : str
        Controls how vertices to be trimmed is determined. Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'full'       Any vertices lying outside the domain are trimmed
            'mixed'      Vertices with at least one neighbor lying inside
                         the domain are kept
            ===========  =====================================================

    """
    from openpnm.topotools import isoutside
    Pdrop = isoutside(network['vert.coords'], shape=shape, thresh=0)
    if mode == 'full':
        network = trim(network=network, pores=np.where(Pdrop)[0])
    elif mode == 'mixed':
        # Find throats connecting internal to external pores
        Ts = np.sum(Pdrop[network['edge.conns']], axis=1) == 1
        # Keep the pores on the ends of these throats
        Pkeep = np.unique(network['edge.conns'][Ts])
        Ps = np.array(list(set(np.where(Pdrop)[0]).difference(set(Pkeep)))).astype(int)
        network = trim(network=network, pores=Ps)
        # Remove throats between these surviving external pores
        Pdrop = isoutside(network['vert.coords'], shape=shape)
        Ts = np.all(Pdrop[network['edge.conns']], axis=1)
        network = trim(network=network, throats=Ts)
        # Lastly label the surviving pores as outside for further processing
        network['vert.outside'] = np.zeros(network['vert.coords'].shape[0], dtype=bool)
        network['vert.outside'][Pdrop] = True
    return network


def _template_sphere_disc(dim, outer_radius, inner_radius):
    r"""
    This private method generates an image array of a sphere/shell-disc/ring.

    It is useful for passing to Cubic networks as a ``template`` to make
    networks with desired shapes.

    Parameters
    ----------
    dim : int
        Network dimension
    outer_radius : int
        Number of the nodes in the outer radius of the network
    inner_radius : int
        Number of the nodes in the inner radius of the network

    Returns
    -------
    im : array_like
        A Numpy array containing 1's to demarcate the desired shape, and 0's
        elsewhere.

    """
    rmax = np.array(outer_radius, ndmin=1)
    rmin = np.array(inner_radius, ndmin=1)
    ind = 2 * rmax - 1
    coord = np.indices((ind * np.ones(dim, dtype=int)))
    coord = coord - (ind - 1)/2
    x = coord[0, :]
    y = coord[1, :]
    if dim == 2:
        img = (x ** 2 + y ** 2) < rmax ** 2
    elif dim == 3:
        z = coord[2, :]
        img = (x ** 2 + y ** 2 + z ** 2) < rmax ** 2
    if rmin[0] != 0:
        if dim == 2:
            img_min = (x ** 2 + y ** 2) > rmin ** 2
        elif dim == 3:
            img_min = (x ** 2 + y ** 2 + z ** 2) > rmin ** 2
        img = img * img_min
    return img


def template_sphere_shell(outer_radius, inner_radius=0, dim=3):
    r"""
    This method generates an image array of a sphere-shell.

    It is useful for passing to Cubic networks as a ``template`` to make
    spherical shaped networks.

    Parameters
    ----------
    outer_radius : int
        Number of nodes in the outer radius of the sphere.
    inner_radius : int
        Number of nodes in the inner radius of the shell.  a value of 0 will
        result in a solid sphere.
    dim : scalar
        Controls the number of dimensions of the result.  3 returns a sphere,
        while 2 returns a disk.

    Returns
    -------
    im : array_like
        A Numpy array containing 1's to demarcate the sphere-shell, and 0's
        elsewhere.

    """
    img = _template_sphere_disc(dim=dim, outer_radius=outer_radius,
                                inner_radius=inner_radius)
    return img


def template_cylinder_annulus(height, outer_radius, inner_radius=0):
    r"""
    This method generates an image array of a disc-ring.

    It is useful for passing to Cubic networks as a ``template`` to make
    circular-shaped 2D networks.

    Parameters
    ----------
    height : int
        The height of the cylinder
    outer_radius : int
        Number of nodes in the outer radius of the cylinder
    inner_radius : int
        Number of the nodes in the inner radius of the annulus.  A value of 0
        will result in a solid cylinder.

    Returns
    -------
    im : array_like
        A Numpy array containing 1's to demarcate the disc-ring, and 0's
        elsewhere.

    """

    img = _template_sphere_disc(dim=2, outer_radius=outer_radius,
                                inner_radius=inner_radius)
    img = np.tile(np.atleast_3d(img), reps=height)
    return img


def reflect_base_points(base_pts, domain_size):
    r"""
    Helper function for relecting a set of points about the faces of a
    given domain.

    Parameters
    ----------
    base_pts : array_like
        The coordinates of the base_pts to be reflected in the coordinate
        system corresponding to the the domain as follows:

        **spherical** : [r, theta, phi]
        **cylindrical** or **circular** : [r, theta, z]
        **rectangular** or **square** : [x, y, z]
    domain_size : list or array
        Controls the size and shape of the domain, as follows:

        **sphere** : If a single value is received, its treated as the radius
        [r] of a sphere centered on [0, 0, 0].

        **cylinder** : If a two-element list is received it's treated as the
        radius and height of a cylinder [r, z] positioned at [0, 0, 0] and
        extending in the positive z-direction.  If the z dimension is 0, a
        disk of radius r is created.

        **rectangle** : If a three element list is received, it's treated
        as the outer corner of rectangle [x, y, z] whose opposite corner lies
        at [0, 0, 0].  If the z dimension is 0, a rectangle of size X-by-Y is
        created.

    Notes
    -----
    The base points can be either [N x 3] or [3 x N].  There transposed internally
    as needed and returned to the original shape.  If N=3 then the transposing is
    skipped so the user needs to ensure the the form of [3 x N].

    """
    domain_size = np.array(domain_size)
    if len(domain_size) == 1:
        r, theta, phi = base_pts
        new_r = 2*domain_size[0] - r
        r = np.hstack([r, new_r])
        theta = np.hstack([theta, theta])
        phi = np.hstack([phi, phi])
        base_pts = np.vstack((r, theta, phi))
    if len(domain_size) == 2:
        r, theta, z = base_pts
        new_r = 2*domain_size[0] - r
        r = np.hstack([r, new_r])
        theta = np.hstack([theta, theta])
        z = np.hstack([z, z])
        if domain_size[1] != 0:  # If not a disk
            r = np.hstack([r, r, r])
            theta = np.hstack([theta, theta, theta])
            z = np.hstack([z, -z, 2*domain_size[1]-z])
        base_pts = np.vstack((r, theta, z))
    elif len(domain_size) == 3:
        Nx, Ny, Nz = domain_size
        # Reflect base points about all 6 faces
        orig_pts = base_pts
        base_pts = np.vstack((base_pts,
                              [-1, 1, 1] * orig_pts + [2.0 * Nx, 0, 0]))
        base_pts = np.vstack((base_pts, [-1, 1, 1] * orig_pts))
        base_pts = np.vstack((base_pts,
                              [1, -1, 1] * orig_pts + [0, 2.0 * Ny, 0]))
        base_pts = np.vstack((base_pts, [1, -1, 1] * orig_pts))
        if domain_size[2] != 0:
            base_pts = np.vstack((base_pts,
                                  [1, 1, -1] * orig_pts + [0, 0, 2.0 * Nz]))
            base_pts = np.vstack((base_pts, [1, 1, -1] * orig_pts))
    return base_pts



def generate_base_points(num_points, domain_size, density_map=None,
                         reflect=True):
    r"""
    Generates a set of base points for passing into the Tessellation-based
    Network classes.  The points can be distributed in spherical, cylindrical,
    or rectilinear patterns, as well as 2D and 3D (disks and squares).

    Parameters
    ----------
    num_points : scalar
        The number of base points that lie within the domain.  Note that the
        actual number of points returned will be larger, with the extra points
        lying outside the domain.
    domain_size : list or array
        Controls the size and shape of the domain, as follows:

        **sphere** : If a single value is received, its treated as the radius
        [r] of a sphere centered on [0, 0, 0].

        **cylinder** : If a two-element list is received it's treated as the
        radius and height of a cylinder [r, z] positioned at [0, 0, 0] and
        extending in the positive z-direction.  If the z dimension is 0, a
        disk of radius r is created.

        **rectangle** : If a three element list is received, it's treated
        as the outer corner of rectangle [x, y, z] whose opposite corner lies
        at [0, 0, 0].  If the z dimension is 0, a rectangle of size X-by-Y is
        created.
    density_map : array, optional
        A an array that contains fractional values (0 < i < 1) indicating the
        liklihood that a point in that region should be kept.  The size of this
        array can be anything, but the shape must match the ``domain_size``;
        that is for a 3D network the shape of the ``density_map`` can be
        [10, 10, 10] or [50, 50, 50], depending on how important the resolution
        of the density distribution is.  For a 2D network the ``density_map``
        should be [10, 10].

        When specifying a custom probabiliy map is it recommended to also set
        values outside the given domain to zero.  If not, then the correct
        shape will still be returned, but with too few points in it.
    reflect : bool
        If ``True``, the the base points are generated as specified, the reflected
        about each face of the domain.  This essentially tricks the
        tessellation functions into creating smoothfaces at the
        boundaries once these excess pores are trimmed.

    Notes
    -----
    The reflection approach tends to create larger pores near the surfaces, so
    it might be necessary to use the ``density_map`` argument to specify a
    slightly higher density of points near the surfaces.

    The ``Voronoi``, ``Delaunay``, ``Gabriel``, and ``DelunayVoronoiDual``
    classes can *techncially* handle base points with spherical or cylindrical
    domains, but the reflection across round surfaces does not create perfect
    Voronoi cells so the surfaces will not be smooth.


    Examples
    --------
    The following generates a spherical array with higher values near the core.
    It uses a distance transform to create a sphere of radius 10, then a
    second distance transform to create larger values in the center away from
    the sphere surface.  These distance values could be further skewed by
    applying a power, with values higher than 1 resulting in higher values in
    the core, and fractional values smoothinging them out a bit.

    >>> import openpnm as op
    >>> import scipy as sp
    >>> import scipy.ndimage as spim
    >>> im = np.ones([21, 21, 21], dtype=int)
    >>> im[10, 10, 10] = 0
    >>> im = spim.distance_transform_edt(im) <= 20  # Create sphere of 1's
    >>> prob = spim.distance_transform_edt(im)
    >>> prob = prob / np.amax(prob)  # Normalize between 0 and 1
    >>> pts = op.topotools.generate_base_points(num_points=50,
    ...                                         domain_size=[1, 1, 1],
    ...                                         density_map=prob)
    >>> net = op.network.DelaunayVoronoiDual(points=pts, shape=[1, 1, 1])

    """

    def _try_points(num_points, prob):
        prob = np.atleast_3d(prob)
        prob = np.array(prob)/np.amax(prob)  # Ensure prob is normalized
        base_pts = []
        N = 0
        while N < num_points:
            pt = np.random.rand(3)  # Generate a point
            # Test whether to keep it or not
            [indx, indy, indz] = np.floor(pt*np.shape(prob)).astype(int)
            if np.random.rand(1) <= prob[indx][indy][indz]:
                base_pts.append(pt)
                N += 1
        base_pts = np.array(base_pts)
        return base_pts

    if len(domain_size) == 1:  # Spherical
        domain_size = np.array(domain_size)
        r = domain_size[0]
        if density_map is None:
            # Make an image of a sphere filled with ones and use _try_points
            density_map = np.ones([41, 41, 41])
            density_map[20, 20, 20] = 0
            density_map = spim.distance_transform_edt(density_map) < 20
        base_pts = _try_points(num_points, density_map)
        # Convert to spherical coordinates
        X, Y, Z = np.array(base_pts - [0.5, 0.5, 0.5]).T
        r = 2*np.sqrt(X**2 + Y**2 + Z**2)*domain_size[0]
        theta = 2*np.arctan(Y/X)
        phi = 2*np.arctan(np.sqrt(X**2 + Y**2)/Z)
        # Trim points outside the domain (from improper prob images)
        inds = r <= domain_size[0]
        [r, theta, phi] = [r[inds], theta[inds], phi[inds]]
        # Reflect base points across perimeter
        if reflect:
            r, theta, phi = reflect_base_points(np.vstack((r, theta, phi)),
                                                domain_size)
        # Convert to Cartesean coordinates
        X, Y, Z = from_spherical(r, theta, phi)
        base_pts = np.vstack([X, Y, Z]).T

    elif len(domain_size) == 2:  # Cylindrical or Disk
        domain_size = np.array(domain_size)
        if density_map is None:
            density_map = np.ones([41, 41, 41])
            density_map[20, 20, :] = 0
            if domain_size[1] == 0:  # Disk
                density_map = density_map[:, :, 0]
            density_map = spim.distance_transform_edt(density_map) < 20
        base_pts = _try_points(num_points, density_map)
        # Convert to cylindrical coordinates
        X, Y, Z = np.array(base_pts - [0.5, 0.5, 0]).T  # Center on z-axis
        r = 2*np.sqrt(X**2 + Y**2)*domain_size[0]
        theta = 2*np.arctan(Y/X)
        z = Z*domain_size[1]
        # Trim points outside the domain (from improper prob images)
        inds = r <= domain_size[0]
        [r, theta, z] = [r[inds], theta[inds], z[inds]]
        inds = ~((z > domain_size[1]) + (z < 0))
        [r, theta, z] = [r[inds], theta[inds], z[inds]]
        if reflect:
            r, theta, z = reflect_base_points(np.vstack([r, theta, z]),
                                              domain_size)
        # Convert to Cartesean coordinates
        X, Y, Z = from_cylindrical(r, theta, z)
        base_pts = np.vstack([X, Y, Z]).T

    elif len(domain_size) == 3:  # Cube or square
        if density_map is None:
            density_map = np.ones([41, 41, 41])
            if domain_size[2] == 0:
                density_map = density_map[:, :, 0]
        base_pts = _try_points(num_points, density_map)
        base_pts = base_pts*domain_size
        if reflect:
            base_pts = reflect_base_points(base_pts, domain_size)

    return base_pts
