import numpy as np
from scipy.spatial import ConvexHull
from transforms3d import _gohlketransforms as tr
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


def PolyArea2D(pts):
    r"""
    returns the area of a 2D polygon given the set of points defining the convex hull
    in correct order
    Example
    ---------
    >>> import OpenPNM.Utilities.vertexops as vo
    >>> tri = np.array([[0,0],[1,2],[2,0]])
    >>> vo.PolyArea2D(tri) == 2.0
    True
    """
    lines = np.hstack([pts, np.roll(pts, -1, axis=0)])
    area = 0.5*abs(sum(x1*y2-x2*y1 for x1, y1, x2, y2 in lines))
    return area


def PolyPerimeter2D(pts):
    r"""
    returns the perimeter of a 2D polygon given the set of points defining the convex
    hull in correct order
    Example
    ---------
    >>> import OpenPNM.Utilities.vertexops as vo
    >>> quad = np.array([[0,0],[0,1],[1,1],[1,0]])
    >>> vo.PolyPerimeter2D(quad) == 4.0
    True
    """
    lines = np.hstack([pts, np.roll(pts, -1, axis=0)])
    perimeter = sum(np.sqrt((x2-x1)**2+(y2-y1)**2) for x1, y1, x2, y2 in lines)
    return perimeter


def PolyWeightedCentroid2D(pts):
    r"""
    returns the centroid of a 2D polygon given the set of points defining the convex
    hull in correct order
    Example
    ---------
    >>> import OpenPNM.Utilities.vertexops as vo
    >>> quad = np.array([[0,0],[0,2],[2,2],[2,0]])
    >>> vo.PolyWeightedCentroid2D(quad) == [1.0,1.0]
    True
    """
    lines = np.hstack([pts, np.roll(pts, -1, axis=0)])
    twice_area = 0.0
    cx = 0.0
    cy = 0.0
    for x1, y1, x2, y2 in lines:
        f = x1*y2 - x2*y1
        cx += (x1+x2)*f
        cy += (y1+y2)*f
        twice_area += f
    A = 0.5*(twice_area)
    Cx = cx/(6*A)
    Cy = cy/(6*A)

    return [Cx, Cy]


def scale(network, scale_factor=[1, 1, 1], preserve_vol=False,
          linear_scaling=[False, False, False]):
    r"""
    A method for scaling the coordinates and vertices to create anisotropic networks
    The original domain volume can be preserved by setting preserve_vol = True

    Example
    ---------
    >>> import OpenPNM
    >>> import OpenPNM.Utilities.vertexops as vo
    >>> import numpy as np
    >>> pn = OpenPNM.Network.Delaunay(num_pores=100, domain_size=[3,2,1])
    >>> pn.add_boundaries()
    >>> B1 = pn.pores("left_boundary")
    >>> B2 = pn.pores("right_boundary")
    >>> Vol = vo.vertex_dimension(pn,B1,B2)
    >>> vo.scale(network=pn,scale_factor=[2,1,1],preserve_vol=True)
    >>> Vol2 = vo.vertex_dimension(pn,B1,B2)
    >>> np.around(Vol-Vol2,5) == 0.0
    True
    >>> vo.scale(network=pn,scale_factor=[2,1,1],preserve_vol=False)
    >>> Vol3 = vo.vertex_dimension(pn,B1,B2)
    >>> np.around(Vol3/Vol,5) == 2.0
    True

    """
    from scipy.special import cbrt
    import scipy as sp
    minmax = np.around(vertex_dimension(network=network,
                                        face1=network.pores(), parm='minmax'), 10)
    scale_factor = np.asarray(scale_factor)
    if preserve_vol is True:
        scale_factor = scale_factor/(cbrt(sp.prod(scale_factor)))

    lin_scale = _linear_scale_factor(network["pore.coords"], minmax,
                                     scale_factor, linear_scaling)

    network["pore.coords"] = network["pore.coords"]*lin_scale
    # Cycle through all vertices of all pores updating vertex values
    for pore in network.pores():
        for i, vert in network['pore.vert_index'][pore].items():
            vert_scale = _linear_scale_factor(vert, minmax,
                                              scale_factor, linear_scaling)
            network["pore.vert_index"][pore][i] = vert*vert_scale
    # Cycle through all vertices of all throats updating vertex values
    for throat in network.throats():
        for i, vert in network['throat.vert_index'][throat].items():

            vert_scale = _linear_scale_factor(vert, minmax,
                                              scale_factor, linear_scaling)
            network["throat.vert_index"][throat][i] = vert*vert_scale
    # Scale the vertices on the voronoi diagram stored on the network
    # These are used for adding boundaries on the Delaunay network class
    vert = network._vor.vertices
    vert_scale = _linear_scale_factor(vert, minmax, scale_factor, linear_scaling)
    network._vor.vertices = vert*vert_scale


def _linear_scale_factor(points=None,
                         minmax=[0, 1, 0, 1, 0, 1],
                         scale_factor=[1, 1, 1],
                         linear_scaling=[False, False, False]):
    r"""
    Work out the linear scale factor of a point or set of points based on the
    domain extent, an absolute scale factor and linear_scaling booleans for each
    coordinate. If all False the absolute scaling is applied equally across the
    domain. If one linear_scaling boolean is True then a linear function is
    applied to scaling the co-ordinates along that axis. If more than one boolean
    is true then a combined linear function is applied
    """
    [xmin, xmax, ymin, ymax, zmin, zmax] = minmax
    max_array = np.array([(xmax-xmin), (ymax-ymin), (zmax-zmin)])
    pos_array = points/max_array
    shape = np.shape(points)
    if len(shape) == 1:
        combined_pos = np.ones([1])
        for i in range(3):
            if linear_scaling[i] is True:
                combined_pos *= pos_array[i]
        lin_scale = (scale_factor - 1)*combined_pos + 1
    else:
        combined_pos = np.ones([shape[0]])
        for i in range(3):
            if linear_scaling[i] is True:
                combined_pos *= pos_array[:, i]
        pos_array = np.vstack((combined_pos, combined_pos, combined_pos)).T
        lin_scale = (scale_factor - 1)*pos_array + 1

    return lin_scale


def vertex_dimension(network, face1=[], face2=[], parm='volume'):
    r"""
    Return the domain extent based on the vertices

    This function is better than using the pore coords as they may be far
    away from the original domain size.  And will alter the effective
    properties which should be based on the original domain sizes. Takes
    one or two sets of pores and works out different geometric properties
    if "length" is specified and two lists are given the planarity is
    determined and the appropriate length (x,y,z) is returned.  It should
    work the same as domain length and area if vertices are not in network
    by using coordinates.

    Example
    ----------
    >>> import OpenPNM
    >>> import OpenPNM.Utilities.vertexops as vo
    >>> pn = OpenPNM.Network.Delaunay(num_pores=100, domain_size=[3, 2, 1])
    >>> pn.add_boundaries()
    >>> B1 = pn.pores('left_boundary')
    >>> B2 = pn.pores('right_boundary')
    >>> vo.vertex_dimension(pn, B1, B2,'volume')
    6.0
    >>> vo.vertex_dimension(pn, B1, B2,'area')
    3.0
    >>> vo.vertex_dimension(pn, B1, B2,'length')
    2.0
    >>> vo.vertex_dimension(pn, B1, B2, 'area_xy')
    6.0
    >>> vo.vertex_dimension(pn, B1, B2, 'area_yz')
    2.0
    >>> vo.vertex_dimension(pn, B1, B2, 'area_xz')
    3.0
    >>> vo.vertex_dimension(pn, B1, B2, 'minmax') == [0.0, 3.0, 0.0, 2.0, 0.0, 1.0]
    True
    """
    pores = np.array([], dtype=int)
    if 0 < len(face1):
        pores = np.hstack((pores, face1))
    if 0 < len(face2):
        pores = np.hstack((pores, face2))

    face1_coords = network['pore.coords'][face1]
    face2_coords = network['pore.coords'][face2]
    face1_planar = np.zeros(3)
    face2_planar = np.zeros(3)
    planar = np.zeros(3)
    for i in range(3):
        if len(np.unique(face1_coords[:, i])) == 1:
            face1_planar[i] = 1
        if len(np.unique(face2_coords[:, i])) == 1:
            face2_planar[i] = 1
    if 0 < len(face1) and 0 < len(face2):
        planar = face1_planar*face2_planar
    elif 0 < len(face1):
        planar = face1_planar
    elif 0 < len(face2):
        planar = face2_planar
    else:
        return 0

    if 'pore.vert_index' in network.props():
        verts = []
        for pore in pores:
            for vert in np.asarray(list(network['pore.vert_index'][pore].values())):
                verts.append(vert)
        verts = np.asarray(verts)
    else:
        verts = network['pore.coords'][pores]

    vx_min = verts[:, 0].min()
    vx_max = verts[:, 0].max()
    vy_min = verts[:, 1].min()
    vy_max = verts[:, 1].max()
    vz_min = verts[:, 2].min()
    vz_max = verts[:, 2].max()
    output = 0
    width = np.around(vx_max-vx_min, 10)
    depth = np.around(vy_max-vy_min, 10)
    height = np.around(vz_max-vz_min, 10)

    if parm == 'volume':
        output = width*depth*height
    elif parm == 'area_xy' or (parm == 'area' and planar[2] == 1):
        output = width*depth
    elif parm == 'area_xz' or (parm == 'area' and planar[1] == 1):
        output = width*height
    elif parm == 'area_yz' or (parm == 'area' and planar[0] == 1):
        output = depth*height
    elif parm == 'length_x' or (parm == 'length' and planar[0] == 1):
        output = width
    elif parm == 'length_y'or (parm == 'length' and planar[1] == 1):
        output = depth
    elif parm == 'length_z'or (parm == 'length' and planar[2] == 1):
        output = height
    elif parm == 'minmax':
        output = [vx_min, vx_max, vy_min, vy_max, vz_min, vz_max]

    return output


def porosity(network):
    r"""
    Return the porosity of the domain - sum of the pore volumes divided by domain
    volume
    """
    domain_vol = vertex_dimension(network, network.pores(), parm='volume')
    try:
        pore_vol = np.sum(network['pore.volume'])
    except KeyError:
        logger.error('Geometries must be assigned first')
        pore_vol = 0
    porosity = np.around(pore_vol/domain_vol, 3)
    return porosity


def pore2centroid(network):
    r"""
    Move the pore coordinate to the centroid of the pore vertices
    """
    for geom_name in network.geometries():
        geometry = network.geometries(geom_name)[0]
        if 'pore.centroid' in geometry.props():
            net_pores, geom_pores = geometry.map_pores(network,
                                                       geometry.pores(),
                                                       True).values()
            for i in range(len(geom_pores)):
                network['pore.coords'][net_pores[i]] = \
                    geometry['pore.centroid'][geom_pores[i]]


def tortuosity(network=None):
    r"""
    Calculate the tortuosity from the angle between throat vectors and principle axes
    """
    conns = network['throat.conns']
    va = network['throat.centroid'] - network['pore.centroid'][conns[:, 0]]
    vb = network['throat.centroid'] - network['pore.centroid'][conns[:, 1]]
    x = [1, 0, 0]
    y = [0, 1, 0]
    z = [0, 0, 1]
    f = 180 / np.pi
    theta_x_a = tr.angle_between_vectors(va, x, directed=False, axis=1)
    theta_x_b = tr.angle_between_vectors(vb, x, directed=False, axis=1)
    theta_x = (np.mean(theta_x_a[~np.isnan(theta_x_a)]) +
               np.mean(theta_x_b[~np.isnan(theta_x_b)])) / 2
    theta_y_a = tr.angle_between_vectors(va, y, directed=False, axis=1)
    theta_y_b = tr.angle_between_vectors(vb, y, directed=False, axis=1)
    theta_y = (np.mean(theta_y_a[~np.isnan(theta_y_a)]) +
               np.mean(theta_y_b[~np.isnan(theta_y_b)])) / 2
    theta_z_a = tr.angle_between_vectors(va, z, directed=False, axis=1)
    theta_z_b = tr.angle_between_vectors(vb, z, directed=False, axis=1)
    theta_z = (np.mean(theta_z_a[~np.isnan(theta_z_a)]) +
               np.mean(theta_z_b[~np.isnan(theta_z_b)])) / 2
    tot_angle = (theta_x+theta_y+theta_z)*f
    if 180 < tot_angle:
        logger.error('Something is wrong: ' + str(tot_angle))

    return 1 / np.cos(np.array([theta_x, theta_y, theta_z]))


def plot_throat(geometry, throats, fig=None):
    r"""
    Print a given throat or list of throats accepted as [1, 2, 3, ..., n]
    Original vertices plus offset vertices are rotated to align with
    the z-axis and then printed in 2D
    e.g vo.print_throat(geom, [34, 65, 99])
    """
    import matplotlib.pyplot as plt
    throat_list = []
    for throat in throats:
        if throat in range(geometry.num_throats()):
            throat_list.append(throat)
        else:
            logger.warn('Throat: ' + str(throat) + ' not part of geometry')
    if len(throat_list) > 0:
        verts = geometry['throat.vertices'][throat_list]
        offsets = geometry['throat.offset_vertices'][throat_list]
        normals = geometry['throat.normal'][throat_list]
        coms = geometry['throat.centroid'][throat_list]
        incentre = geometry['throat.incentre'][throat_list]
        inradius = 0.5*geometry['throat.indiameter'][throat_list]
        row_col = np.ceil(np.sqrt(len(throat_list)))
        for i in range(len(throat_list)):
            if fig is None:
                fig = plt.figure()
            ax = fig.add_subplot(row_col, row_col, i+1)
            vert_2D = rotate_and_chop(verts[i], normals[i], [0, 0, 1])
            hull = ConvexHull(vert_2D, qhull_options='QJ Pp')
            for simplex in hull.simplices:
                plt.plot(vert_2D[simplex, 0], vert_2D[simplex, 1], 'k-', linewidth=2)
            plt.scatter(vert_2D[:, 0], vert_2D[:, 1])
            offset_2D = rotate_and_chop(offsets[i], normals[i], [0, 0, 1])
            offset_hull = ConvexHull(offset_2D, qhull_options='QJ Pp')
            for simplex in offset_hull.simplices:
                plt.plot(offset_2D[simplex, 0], offset_2D[simplex, 1],
                         'g-', linewidth=2)
            plt.scatter(offset_2D[:, 0], offset_2D[:, 1])
            # Make sure the plot looks nice by finding the greatest range of points
            # and setting the plot to look square
            xmax = vert_2D[:, 0].max()
            xmin = vert_2D[:, 0].min()
            ymax = vert_2D[:, 1].max()
            ymin = vert_2D[:, 1].min()
            x_range = xmax - xmin
            y_range = ymax - ymin
            if (x_range > y_range):
                my_range = x_range
            else:
                my_range = y_range
            lower_bound_x = xmin - my_range*0.5
            upper_bound_x = xmin + my_range*1.5
            lower_bound_y = ymin - my_range*0.5
            upper_bound_y = ymin + my_range*1.5
            plt.axis((lower_bound_x, upper_bound_x, lower_bound_y, upper_bound_y))
            plt.grid(b=True, which='major', color='b', linestyle='-')
            centroid = rotate_and_chop(coms[i], normals[i], [0, 0, 1])
            incent = rotate_and_chop(incentre[i], normals[i], [0, 0, 1])
            plt.scatter(centroid[0][0], centroid[0][1])
            # Plot incircle
            t = np.linspace(0, 2*np.pi, 200)
            u = inradius[i]*np.cos(t)+incent[0][0]
            v = inradius[i]*np.sin(t)+incent[0][1]
            plt.plot(u, v, 'r-')
            ax.ticklabel_format(style='sci', scilimits=(0, 0))

    else:
        logger.error("Please provide throat indices")
    return fig


def plot_pore(geometry, pores, fig=None, axis_bounds=None, include_points=False):
    r"""
    Print all throats around a given pore or list of pores accepted
    as [1, 2, 3, ..., n]
    e.g vo.print_pore(geom, [34, 65, 99])
    Original vertices plus offset vertices used to create faces and
    then printed in 3D
    To print all pores (n)
    pore_range = np.arange(0,n-1,1)
    vo.print_pore(geom, pore_range)
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    if len(pores) > 0:
        net_pores = geometry.map_pores(geometry._net, pores)
        centroids = geometry['pore.centroid'][pores]
        coords = geometry._net['pore.coords'][net_pores]
        net_throats = geometry._net.find_neighbor_throats(pores=net_pores)
        throats = geometry._net.map_throats(geometry,
                                            net_throats,
                                            return_mapping=True)['target']
        tcentroids = geometry["throat.centroid"][throats]
        # Can't create volume from one throat
        if 1 <= len(throats):
            verts = geometry['throat.vertices'][throats]
            normals = geometry['throat.normal'][throats]
            # Get verts in hull order
            ordered_verts = []
            for i in range(len(verts)):
                vert_2D = rotate_and_chop(verts[i], normals[i], [0, 0, 1])
                hull = ConvexHull(vert_2D, qhull_options='QJ Pp')
                ordered_verts.append(verts[i][hull.vertices])
            offsets = geometry['throat.offset_vertices'][throats]
            ordered_offs = []
            for i in range(len(offsets)):
                offs_2D = rotate_and_chop(offsets[i], normals[i], [0, 0, 1])
                offs_hull = ConvexHull(offs_2D, qhull_options='QJ Pp')
                ordered_offs.append(offsets[i][offs_hull.vertices])
            # Get domain extents for setting axis
            if axis_bounds is None:
                [xmin, xmax, ymin, ymax, zmin, zmax] = \
                    vertex_dimension(geometry._net, pores, parm='minmax')
            else:
                [xmin, xmax, ymin, ymax, zmin, zmax] = axis_bounds
            if fig is None:
                fig = plt.figure()
            ax = fig.gca(projection='3d')
            outer_items = Poly3DCollection(ordered_verts, linewidths=1,
                                           alpha=0.2, zsort='min')
            outer_face_colours = [(1, 0, 0, 0.01)]
            outer_items.set_facecolor(outer_face_colours)
            ax.add_collection(outer_items)
            inner_items = Poly3DCollection(ordered_offs, linewidths=1,
                                           alpha=0.2, zsort='min')
            inner_face_colours = [(0, 0, 1, 0.01)]
            inner_items.set_facecolor(inner_face_colours)
            ax.add_collection(inner_items)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_zlim(zmin, zmax)
            if include_points:
                ax.scatter(centroids[:, 0], centroids[:, 1], centroids[:, 2], c='y')
                ax.scatter(tcentroids[:, 0], tcentroids[:, 1], tcentroids[:, 2],
                           c='r')
                ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c='b')
            ax.ticklabel_format(style='sci', scilimits=(0, 0))

        else:
            plot_throat(geometry, throats, fig)
    else:
        logger.error('Please provide pore indices')
    return fig


def rotate_and_chop(verts, normal, axis=[0, 0, 1]):
    r"""
    Method to rotate a set of vertices (or coords) to align with an axis
    points must be coplanar and normal must be given
    Chops axis coord to give vertices back in 2D
    Used to prepare verts for printing or calculating convex hull in order to arrange
    them in hull order for calculations and printing
    """
    xaxis = [1, 0, 0]
    yaxis = [0, 1, 0]
    zaxis = [0, 0, 1]
    angle = tr.angle_between_vectors(normal, axis)
    if angle == 0.0 or angle == np.pi:
        # We are already aligned
        facet = verts
    else:
        M = tr.rotation_matrix(tr.angle_between_vectors(normal, axis),
                               tr.vector_product(normal, axis))
        try:
            facet = np.dot(verts, M[:3, :3].T)
        except ValueError:
            pass

    try:
        x = facet[:, 0]
        y = facet[:, 1]
        z = facet[:, 2]
    except IndexError:
        x = facet[0]
        y = facet[1]
        z = facet[2]
    # Work out span of points and set axes scales to cover this and be
    # equal in both dimensions
    if axis == xaxis:
        output = np.column_stack((y, z))
    elif axis == yaxis:
        output = np.column_stack((x, z))
    elif axis == zaxis:
        output = np.column_stack((x, y))
    else:
        output = facet

    return output
