"""
Tools
-----

"""
import numpy as np
from openpnm._skgraph import tools


__all__ = [
    'parse_points',
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
    # Deal with input arguments
    if isinstance(points, int):
        points = generate_base_points(num_points=points,
                                      domain_size=shape,
                                      reflect=reflect)
    else:
        # Should we check to ensure that points are reflected?
        points = np.array(points)
    return points


def add_all_label(network):
    r"""
    Add 'node.all' and 'edge.all' to network dictionary

    Parameters
    ----------
    network : dict
        The network dictionary, containing 'node.coords' and 'edge.conns'

    Returns
    -------
    network : dict
        The supplied dictionary with the 'all' labels added

    Notes
    -----
    This function is helpful for working with OpenPNM

    """
    node_prefix = tools.get_node_prefix(network)
    edge_prefix = tools.get_edge_prefix(network)
    coords = network[node_prefix+'.coords']
    conns = network[edge_prefix+'.conns']
    network['pore.all'] = np.ones(coords.shape[0], dtype=bool)
    network['throat.all'] = np.ones(conns.shape[0], dtype=bool)
    return network


def label_faces_cubic(network, rtol=0.0):
    r"""
    Label the nodes sitting on the faces of the domain assuming the domain
    is cubic

    Parameters
    ----------
    network : dict
        The network dictionary contain 'node.coords'
    rtol : float
        Controls how closely a node must be to a face to be counted. It is
        computed relative to the fraction of domain size, as:
        ``hi_label = abs(1 - x[i]/x.max()) < rtol`` and
        ``lo_label = abs(x[i]/x.max()) < rtol``

    Returns
    -------
    network : dict
        The network dictionary with the face labels added

    """
    node_prefix = tools.get_node_prefix(network)
    coords = network[node_prefix+'.coords']
    dims = tools.dimensionality(network)
    coords = np.around(coords, decimals=10)
    min_labels = ['left', 'front', 'bottom']
    max_labels = ['right', 'back', 'top']
    min_coords = np.amin(coords, axis=0)
    max_coords = np.amax(coords, axis=0)
    for ax in np.where(dims)[0]:
        network[node_prefix + '.' + min_labels[ax]] = \
            abs((coords[:, ax]-min_coords[ax])/max_coords[ax]) <= rtol
        network[node_prefix + '.' + max_labels[ax]] = \
            abs(1-coords[:, ax]/max_coords[ax]) <= rtol
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


def template_sphere_shell(r_outer, r_inner=0):
    r"""
    This method generates an image array of a sphere-shell.

    It is useful for passing to Cubic networks as a ``template`` to make
    spherical shaped networks.

    Parameters
    ----------
    r_outer : int
        Number of nodes in the outer radius of the sphere
    r_inner : int, optional
        Number of nodes in the inner radius of the shell.  A value of 0 will
        result in a solid sphere.

    Returns
    -------
    im : array_like
        A Numpy array containing 1's to demarcate the sphere-shell, and 0's
        elsewhere.

    """
    img = _template_sphere_disc(dim=3, outer_radius=r_outer,
                                inner_radius=r_inner)
    return img


def template_cylinder_annulus(z, r_outer, r_inner=0):
    r"""
    This method generates an image array of a disc-ring.

    It is useful for passing to Cubic networks as a ``template`` to make
    circular-shaped 2D networks.

    Parameters
    ----------
    z : int
        The height of the cylinder. A value of 0 will result in a circle
    r_outer : int
        Number of nodes in the outer radius of the cylinder
    r_inner : int, optional
        Number of the nodes in the inner radius of the annulus.  A value of 0
        will result in a solid cylinder.

    Returns
    -------
    im : array_like
        A Numpy array containing 1's to demarcate the disc-ring, and 0's
        elsewhere.

    """

    img = _template_sphere_disc(dim=2, outer_radius=r_outer,
                                inner_radius=r_inner)
    if z > 0:
        img = np.tile(np.atleast_3d(img), reps=z)
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
                   points in the range [-r:r, -r:r, 0:z]
        [r, 0]     A 2D circular domain of radius r, with points in the range
                   [-r:r, -r:r, 0:0]
        [r]        A 3D spherical domain of radius r, with points in the range
                   [-r:r, -r:r, -r:r]
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
        shape = np.array([2*domain_size[0], 2*domain_size[0], 2*domain_size[0]])
        base_pts = np.random.rand(num_points*9, 3)  # Generate more than needed
        # Convert to spherical coordinates
        X, Y, Z = ((np.array(base_pts - [0.5, 0.5, 0.5]))*shape).T
        R, Q, P = tools.cart2sph(X, Y, Z)
        # Keep points outside the domain
        inds = R <= domain_size[0]
        R, Q, P = R[inds], Q[inds], P[inds]
        # Trim internal points to give requested final number
        R, Q, P = R[:num_points], Q[:num_points], P[:num_points]
        # Reflect base points across perimeter
        if reflect:
            R, Q, P = reflect_base_points(np.vstack((R, Q, P)), domain_size)
        # Convert to Cartesean coordinates
        base_pts = np.vstack(tools.sph2cart(R, Q, P)).T

    elif len(domain_size) == 2:  # Cylindrical or Disk
        shape = np.array([2*domain_size[0], 2*domain_size[0], domain_size[1]])
        base_pts = np.random.rand(num_points*9, 3)  # Generate more than needed
        # Convert to cylindrical coordinates
        X, Y, Z = ((np.array(base_pts - [0.5, 0.5, 0]))*shape).T
        R, Q, Z = tools.cart2cyl(X, Y, Z)
        # Trim points outside the domain
        inds = R <= domain_size[0]
        R, Q, Z = R[inds], Q[inds], Z[inds]
        # Reduce to requested number of points
        R, Q, Z = R[:num_points], Q[:num_points], Z[:num_points]
        if reflect:
            R, Q, Z = reflect_base_points(np.vstack((R, Q, Z)), domain_size)
        # Convert to Cartesean coordinates
        base_pts = np.vstack(tools.cyl2cart(R, Q, Z)).T

    elif len(domain_size) == 3:  # Cube or square
        base_pts = np.random.rand(num_points, 3)*domain_size
        if reflect:
            base_pts = reflect_base_points(base_pts, domain_size)

    return base_pts
