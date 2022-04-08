"""
Tools
-----

"""
import numpy as np
from openpnm._skgraph import settings


__all__ = [
    'parse_points',
    'get_spacing',
    'get_shape',
    'add_all_label',
    'label_surface_nodes',
    'label_faces',
    'template_sphere_shell',
    'template_cylinder_annulus',
]


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


def get_spacing(coords, conns):
    r"""
    Determine spacing of a cubic network

    Parameters
    ----------
    coords : ndarray
        The x,y,z, location of each node
    conns : ndarray
        The pore1-pore2 connections for each throat

    Returns
    -------
    spacing : ndarray
        An array containing the spacing between nodes in each direction

    Notes
    -----
    This function only works on simple cubic networks with no boundary
    vertices. If a unique spacing cannot be found in each direction,
    and/or the edges are not all oriented perpendicularly, exceptions
    will be raised.

    """
    from openpnm.topotools import dimensionality
    # Find Network spacing
    C12 = coords[conns]
    mag = np.linalg.norm(np.diff(C12, axis=1), axis=2)
    unit_vec = np.around(np.squeeze(np.diff(C12, axis=1)) / mag, decimals=14)
    spacing = [0, 0, 0]
    dims = dimensionality(coords=coords)
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


def get_shape(coords, conns):
    L = np.ptp(coords, axis=0)
    mask = L.astype(bool)
    S = get_spacing(coords, conns)
    shape = np.array([1, 1, 1], int)
    shape[mask] = L[mask] / S[mask] + 1
    return shape


def add_all_label(coords, conns):
    d = {}
    d['pore.all'] = np.ones(coords.shape[0], dtype=bool)
    d['throat.all'] = np.ones(conns.shape[0], dtype=bool)
    return d


def label_surface_nodes(coords):
    r"""
    """
    from openpnm._skimage.tools import dimensionality

    node_prefix = settings.node_prefix

    dims = dimensionality(coords)
    hits = np.zeros_like(coords.shape[0], dtype=bool)
    mn = np.amin(coords, axis=0)
    mx = np.amax(coords, axis=0)
    for ax in np.where(dims)[0]:
        if dims[ax]:
            hits += coords[:, ax] <= mn[ax]
            hits += coords[:, ax] >= mx[ax]
    d = {node_prefix + '.surface': hits}
    return d


def label_faces(coords, threshold=0.05):
    r"""
    Label the vertices sitting on the faces of the domain in accordance with
    the conventions used for cubic etc.
    """
    from openpnm._skimage.tools import dimensionality

    node_prefix = settings.node_prefix

    dims = dimensionality(coords)
    coords = np.around(coords, decimals=10)
    min_labels = ['front', 'left', 'bottom']
    max_labels = ['back', 'right', 'top']
    min_coords = np.amin(coords, axis=0)
    max_coords = np.amax(coords, axis=0)
    d = {}
    for ax in np.where(dims)[0]:
        d[node_prefix + '.' + min_labels[ax]] = \
            coords[:, ax] <= threshold*min_coords[ax]
        d[node_prefix + '.' + max_labels[ax]] = \
            coords[:, ax] >= (1-threshold)*max_coords[ax]
    return d


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


def reflect_base_points(points, domain_size):
    r"""
    Helper function for relecting a set of points about the faces of a
    given domain.

    Parameters
    ----------
    points : array_like
        The coordinates of the base_pts to be reflected in the coordinate
        system corresponding to the the domain as follows:

        **spherical** : [r, theta, phi]
        **cylindrical** or **circular** : [r, theta, z]
        **rectangular** or **square** : [x, y, z]
    domain_size : list or array
        Controls the size and shape of the domain, as follows:

        ========== ============================================================
        shape      result
        ========== ============================================================
        [x, y, z]  A 3D cubic domain of dimension x, y and z
        [x, y, 0]  A 2D square domain of size x by y
        [r, z]     A 3D cylindrical domain of radius r and height z
        [r, 0]     A 2D circular domain of radius r
        [r]        A 3D spherical domain of radius r
        ========== ============================================================

    """
    domain_size = np.array(domain_size)
    if len(domain_size) == 1:
        r, theta, phi = points
        new_r = 2*domain_size[0] - r
        r = np.hstack([r, new_r])
        theta = np.hstack([theta, theta])
        phi = np.hstack([phi, phi])
        points = np.vstack((r, theta, phi))
    if len(domain_size) == 2:
        r, theta, z = points
        new_r = 2*domain_size[0] - r
        r = np.hstack([r, new_r])
        theta = np.hstack([theta, theta])
        z = np.hstack([z, z])
        if domain_size[1] != 0:  # If not a disk
            r = np.hstack([r, r, r])
            theta = np.hstack([theta, theta, theta])
            z = np.hstack([z, -z, 2*domain_size[1]-z])
        points = np.vstack((r, theta, z))
    elif len(domain_size) == 3:
        Nx, Ny, Nz = domain_size
        # Reflect base points about all 6 faces
        orig_pts = points
        points = np.vstack((points,
                            [-1, 1, 1] * orig_pts + [2.0 * Nx, 0, 0]))
        points = np.vstack((points, [-1, 1, 1] * orig_pts))
        points = np.vstack((points,
                            [1, -1, 1] * orig_pts + [0, 2.0 * Ny, 0]))
        points = np.vstack((points, [1, -1, 1] * orig_pts))
        if domain_size[2] != 0:
            points = np.vstack((points,
                               [1, 1, -1] * orig_pts + [0, 0, 2.0 * Nz]))
            points = np.vstack((points, [1, 1, -1] * orig_pts))
    return points


def generate_base_points(num_points, domain_size, reflect=True):
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

        ========== ============================================================
        shape      result
        ========== ============================================================
        [x, y, z]  A 3D cubic domain of dimension x, y and z
        [x, y, 0]  A 2D square domain of size x by y
        [r, z]     A 3D cylindrical domain of radius r and height z
        [r, 0]     A 2D circular domain of radius r
        [r]        A 3D spherical domain of radius r
        ========== ============================================================

    reflect : bool
        If ``True``, the base points are generated as specified then reflected
        about each face of the domain.  This essentially tricks the
        tessellation functions into creating smooth faces at the
        boundaries once these excess points are trimmed.

    Notes
    -----
    The ``Voronoi``, ``Delaunay``, ``Gabriel``, and ``DelunayVoronoiDual``
    classes can *techncially* handle base points with spherical or cylindrical
    domains, but the reflection across round surfaces does not create perfect
    Voronoi cells so the surfaces will not be smooth.


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
        base_pts = _try_points(num_points, density_map)
        base_pts = base_pts*domain_size
        if reflect:
            base_pts = reflect_base_points(base_pts, domain_size)

    return base_pts
