"""
Tools
-----

"""
import numpy as np
from openpnm._skgraph import settings
from openpnm._skgraph import tools


__all__ = [
    'parse_points',
    'get_spacing',
    'get_shape',
    'add_all_label',
    'label_faces_cubic',
    'template_sphere_shell',
    'template_cylinder_annulus',
    'generate_base_points',
    'reflect_base_points',
]

__notyet__ = [
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


def label_faces_cubic(coords, threshold=0.0):
    r"""
    Label the nodes sitting on the faces of the domain assuming the domain
    is cubic

    Parameters
    ----------
    coords : ndarray
        The x,y,z coordinates of each node
    threshold : float
        Controls how closely a node must be to a face to be counted. It is
        defined as ``hit = x <= (1+threshold) * x_min`` and
        ``hit = y >= (1-threshold) * y_max``.  The default is 0, which means
        nodes must be exactly on a face.
    """
    node_prefix = settings.node_prefix

    dims = tools.dimensionality(coords)
    coords = np.around(coords, decimals=10)
    min_labels = ['front', 'left', 'bottom']
    max_labels = ['back', 'right', 'top']
    min_coords = np.amin(coords, axis=0)
    max_coords = np.amax(coords, axis=0)
    d = {}
    for ax in np.where(dims)[0]:
        d[node_prefix + '.' + min_labels[ax]] = \
            coords[:, ax] <= (1+threshold)*min_coords[ax]
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
    Relects a set of points about the faces of a given domain

    Parameters
    ----------
    points : ndarray
        The coordinates of the points to be reflected.  The points should be
        in the coordinate system corresponding to the the domain.
    domain_size : list or array
        Controls the size and shape of the domain, as follows:

        ========== ============================================================
        shape      result
        ========== ============================================================
        [x, y, z]  Points will be reflected about all 6 faces
        [x, y, 0]  Points will be relfected about the x and y faces
        [r, z]     Points will be reflected above and below the end faces and
                   across the perimeter
        [r, 0]     Points will be reflected across the perimeter
        [r]        Points will be reflected across the outer surface
        ========== ============================================================

    Returns
    -------
    points : ndarray
        The coordinates of the original points plus the new reflected points
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
    Generates a set of randomly distributed points in rectilinear coordinates
    for use in spatial tessellations

    The points can be distributed in spherical, cylindrical, or rectilinear
    domains, as well as 2D and 3D (disks and squares)

    Parameters
    ----------
    num_points : scalar
        The number of base points that lie within the domain
    domain_size : list or array
        Controls the size and shape of the domain, as follows:

        ========== ============================================================
        shape      result
        ========== ============================================================
        [x, y, z]  A 3D cubic domain of dimension x, y and z, with points in
                   the range [0:x, 0:y, 0:z]
        [x, y, 0]  A 2D square domain of size x by y, with points in the range
                   [0:x, 0:y, 0:0]
        [r, z]     A 3D cylindrical domain of radius r and height z, with
                   points in the range [-r/2:r/2, -r/2:r/2, 0:z]
        [r, 0]     A 2D circular domain of radius r, with points in the range
                   [-r/2:r/2, -r/2:r/2, 0:0]
        [r]        A 3D spherical domain of radius r, with points in the range
                   [-r/2:r/2, -r/2:r/2, -r/2:r/2]
        ========== ============================================================

    reflect : bool
        If ``True``, the base points are generated as specified then reflected
        about each face of the domain.  This essentially tricks the
        tessellation functions into creating smooth faces at the boundaries
        once these excess points are trimmed. Note that the surface is not
        perfectly smooth for the curved faces.

    Notes
    -----
    To convert between coordinate systems see ``cart2sph`` and ``cart2cyl``

    """
    if len(domain_size) == 1:  # Spherical
        domain_size = np.array(domain_size)
        base_pts = np.random.rand(num_points*9, 3)  # Generate more than needed
        # Convert to spherical coordinates
        X, Y, Z = (np.array(base_pts - [0.5, 0.5, 0.5]))*domain_size
        R, Q, P = tools.cart2sph(X, Y, Z).T
        # Keep points outside the domain
        inds = R <= domain_size[0]/2
        R, Q, P = R[inds], Q[inds], P[inds]
        # Trim internal points to give requested final number
        R, Q, P = R[:num_points], Q[:num_points], P[:num_points]
        # Reflect base points across perimeter
        if reflect:
            R, Q, P = reflect_base_points(np.vstack((R, Q, P)).T, domain_size)
        # Convert to Cartesean coordinates
        base_pts = tools.sph2cart(R, Q, P)

    elif len(domain_size) == 2:  # Cylindrical or Disk
        domain_size = np.array([domain_size[0], domain_size[0], domain_size[1]])
        base_pts = np.random.rand(num_points*9, 3)  # Generate more than needed
        # Convert to cylindrical coordinates
        X, Y, Z = ((np.array(base_pts - [0.5, 0.5, 0]))*domain_size).T
        R, Q, Z = tools.cart2cyl(X, Y, Z).T
        # Trim points outside the domain (from improper prob images)
        inds = R <= domain_size[0]/2
        R, Q, Z = R[inds], Q[inds], Q[inds]
        R, Q, Z = R[:num_points], Q[:num_points], Z[:num_points]
        if reflect:
            R, Q, Z = reflect_base_points(np.vstack((R, Q, Z)).T, domain_size)
        # Convert to Cartesean coordinates
        base_pts = tools.cyl2cart(R, Q, Z)

    elif len(domain_size) == 3:  # Cube or square
        base_pts = np.random.rand(num_points, 3)*domain_size
        if reflect:
            base_pts = reflect_base_points(base_pts, domain_size)

    return base_pts
