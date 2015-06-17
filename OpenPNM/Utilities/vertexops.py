import numpy as np
from scipy.spatial import ConvexHull
from OpenPNM.Utilities import transformations as tr
from OpenPNM.Utilities import misc as misc
from math import atan2


def get_throat_geom(verts, normal, fibre_rad):
    r"""
    For one set of vertices defining a throat return the key properties
    This is the main loop for calling other sub-routines.
    General Method:
        For each connection or throat defined by the shared vertices
        Rotate the vertices to align with the xy-plane and get rid of z-coordinate
        Compute the convex hull of the 2D points giving a set of simplices which
        define neighbouring vertices in a clockwise fashion.
        For each triplet calculate the offset position given the fibre radius
        Check for overlapping vertices and ones that lie outside the original hull -
        recalculate position to offset from or ignore if all overlapping.
        Calculate Area and Perimeter if successfully generated offset vertices to
        replicate eroded throat.
        Translate back into 3D
    Any Errors encountered result in the throat area being zero and no vertices being
    passed back.
    These Errors are not coding mistakes but failures to obtain an eroded facet with
    non-zero area:
    Error 1: Less than 3 vertices in the convex hull - Should never happen
             (unless 2 points are incredibly close together)
    Error 2: The largest span of points is less than twice the fibre radius
             (i.e. throat will definitley be occluded)
    Error 3: All the offset vertices overlap with at least one other vertex -
             Throat fully occluded
    Error 4: Not enough offset vertices to continue - Throat fully occluded
    Error 5: An offset vertex is outside the original set of points - Throat fully
             occluded
    """
    z_axis = [0, 0, 1]
    throat_area = 0.0
    throat_perimeter = 0.0
    throat_COM = np.zeros([1, 3])
    output_offset = []
    Error = 0
    # For boundaries some facets will already be aligned with the axis -
    # if this is the case a rotation is unnecessary and could also cause problems
    angle = tr.angle_between_vectors(normal, z_axis)
    if angle == 0.0 or angle == np.pi:
        # We are already aligned
        rotate_input = False
        facet = verts
    else:
        rotate_input = True
        M = tr.rotation_matrix(tr.angle_between_vectors(normal, z_axis),
                               tr.vector_product(normal, z_axis))
        facet = np.dot(verts, M[:3, :3].T)
    x = facet[:, 0]
    y = facet[:, 1]
    z = facet[:, 2]
    # Work out span of points and set axes scales to cover this and be equal
    # in both dimensions
    x_range = x.max() - x.min()
    y_range = y.max() - y.min()
    if (x_range > y_range):
        my_range = x_range
    else:
        my_range = y_range
    if np.around(z.std(), 3) != 0.000:
        print('Rotation failed')
    facet_coords_2D = np.column_stack((x, y))
    hull = ConvexHull(facet_coords_2D, qhull_options='QJ Pp')
    verts_2D = facet_coords_2D[hull.vertices]
    offset = outer_offset(verts_2D, fibre_rad)
    # At this point we may have overlapping areas for which we need to offset
    # from a new point
    overlap_array, sweep_radius, line_points = set_overlap(verts_2D, offset)
    temp_vert_list = []
    if len(verts_2D) < 3:
        # Error: Fused Too Many Verts
        Error = 1
    elif(my_range < fibre_rad*2):
        # Error: Facet Too small to Erode
        Error = 2
    else:
        if all_overlap(overlap_array) is True:
            # If all overlaps then throat is fully occluded
            # Error: Throat fully occluded
            Error = 3
        else:
            # If one or two sets of overlaps exist and at least one vertex is not
            # overlapped then we need to do a bit more work. Do some linalg to find
            # a new point to offset from saving un-overlapped verts and newly
            # created verts in a temporary list
            count = 0
            temp_verts = verts_2D
            while True:
                temp_vert_list = []
                for i in range(np.shape(line_points)[0]):
                    if np.sum(overlap_array[i]) == 0.0:
                        temp_vert_list.append(temp_verts[i])
                    else:
                        my_lines = []
                        for j in range(np.shape(line_points)[0]):

                            if overlap_array[i][j] == 1 and overlap_array[j][i] == 1:
                                list_a = line_points[i][j]
                                list_b = line_points[j][i]
                                my_lines = symmetric_difference(list_a, list_b)

                        my_lines = np.asarray(my_lines)

                        if len(my_lines) == 2:
                            try:
                                quad_points = temp_verts[my_lines]
                                my_new_point = new_point(quad_points)
                                temp_vert_list.append(my_new_point)
                            except IndexError:
                                print('IndexError: ' + str(my_lines))
                            except TypeError:
                                print('TypeError: ' + str(my_lines))

                temp_verts = np.asarray(misc.unique_list(temp_vert_list))
                offset = outer_offset(temp_verts, fibre_rad)
                overlap_array, sweep_radius, line_points = \
                    set_overlap(temp_verts, offset)

                if overlap_array.any() is False:
                    break
                elif all_overlap(overlap_array) is True:
                    Error = 3
                    break
                elif len(temp_verts) < 3:
                    Error = 4
                    break
                else:
                    count += 1
                    temp_verts = np.asarray(fuse_verts(verts=temp_verts,
                                                       percentage=0.05*count))
                    offset = outer_offset(temp_verts, fibre_rad)
                    overlap_array, sweep_radius, line_points = \
                        set_overlap(temp_verts, offset)
                    # Continue looping until one of the above conditions is true
                    # or counter reaches 10
                if count >= 10:
                    break

    if len(offset) >= 3 and Error == 0:
        # Now also check whether any of the offset points lie outside the
        # original convex hull
        original_area = np.around(PolyArea2D(verts_2D), 10)
        all_points = np.concatenate((verts_2D, offset), axis=0)
        try:
            total_hull = ConvexHull(all_points, qhull_options='QJ Pp')
            total_area = np.around(PolyArea2D(all_points[total_hull.vertices]), 10)
        except np.spatial.qhull.QhullError:
            print(all_points)
            total_area = 999
            Error = 5
        offset_hull = ConvexHull(offset, qhull_options='QJ Pp')
        offset_verts_2D = offset[offset_hull.vertices]
        if total_area > original_area:
            if Error != 5:
                Error = 6
        else:
            throat_area = PolyArea2D(offset_verts_2D)
            throat_perimeter = PolyPerimeter2D(offset_verts_2D)
            throat_COM_2D = PolyWeightedCentroid2D(offset_verts_2D)
            throat_COM = np.hstack((throat_COM_2D, z[0]))
        # Make 3D again in rotated plane
        offset_verts_3D = np.column_stack((offset_verts_2D,
                                           z[0:len(offset_verts_2D)]))
        # Get matrix to un-rotate the co-ordinates back to the original orientation
        # if we rotated in the first place
        if (rotate_input):
            M1 = tr.inverse_matrix(M)
            # Unrotate the offset coordinates
            output_offset = np.dot(offset_verts_3D, M1[:3, :3].T)
            throat_COM = np.dot(throat_COM, M1[:3, :3].T)
        else:
            output_offset = offset_verts_3D

    return throat_area, throat_perimeter, output_offset, throat_COM, Error


def outer_offset(verts, fibre_rad):
    r"""
    Routine to loop through all verts and calculate offset position based on
    neighbours either side. Verts must be in hull order
    """
    offset = []
    for i, vert in enumerate(verts):
        # Collect three adjacent points and compute the offset of the first
        triplet = (vert, np.roll(verts, -1, axis=0)[i], np.roll(verts, 1, axis=0)[i])
        offset.append(offset_vertex(triplet, fibre_rad))
    offset = np.asarray(offset)

    return offset


def offset_vertex(points, rad=0.01):
    # We are passed in a set of 3 points forming vertices of two adjoining simplexes
    # of the convex hull of a voronoi facet. We need to offset the vertices normal to
    # the fibre direction (or adjoining vectors) by the fibre radius. This is
    # achieved by finding the half angle between the two adjoining vectors and a
    # direction. Mid-point must be the first in the array
    p0 = np.array(points[0])
    p1 = np.array(points[1])
    p2 = np.array(points[2])
    # Now make the midpoint the origin
    vector1 = p1-p0
    vector2 = p2-p0

    # Find what quadrant the vector is pointing in - atan2 function takes account of
    # signs. 0 means aligned with x-axis, pi is aligned with -xaxis, positive numbers
    # are positive y and negative numbers are negative y. The angle between the
    # vectors should always be within 180 degrees of one another in a convex hull.

    q1 = atan2(vector1[1], vector1[0])
    q2 = atan2(vector2[1], vector2[0])
    alpha = 0.5*tr.angle_between_vectors(vector1, vector2)

    # We always want to offset from the first vertex we get to - going anti-clockwise
    # from the x-axis Check if both vectors point up or both point down - if so the
    # first one we get to will have smaller q value
    if q1*q2 >= 0.0:
        if q1 < q2:
            theta = q1
        else:
            theta = q2
    else:
        # If vector 1 is more rotated away from positive xaxis than vector 2 is
        # rotated away from negative xaxis - use it, and vice-versa
        if (abs(q1) + abs(q2)) > np.pi:
            # Vectors are pointing negative x so take whichever has positive q-value
            # like a pacman facing left
            if q1 >= 0:
                theta = q1
            else:
                theta = q2
        else:
            # Vectors are pointing positive x so take whichever is negative
            if q1 <= 0:
                theta = q1
            else:
                theta = q2

    # This may cause problems in terms of which way to offset.
    if alpha == 0:
        x = 0
        y = 0
    else:
        x = rad*np.cos(alpha+theta)/np.sin(alpha)
        y = rad*np.sin(alpha+theta)/np.sin(alpha)

    # Add the midpoint back in
    output = [x+p0[0], y+p0[1]]

    return output


def dist2(p1, p2):
    r"""
    Pythagoras to compute the square of the distance between two points (in 2D)
    """
    return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2


def fuse(points, d):
    r"""
    Fuse points together wihin a certain range
    """
    ret = []
    d2 = d * d
    n = len(points)
    taken = [False] * n
    for i in range(n):
        if not taken[i]:
            count = 1
            point = [points[i][0], points[i][1]]
            taken[i] = True
            for j in range(i+1, n):
                if dist2(points[i], points[j]) < d2:
                    point[0] += points[j][0]
                    point[1] += points[j][1]
                    count += 1
                    taken[j] = True
            point[0] /= count
            point[1] /= count
            ret.append((point[0], point[1]))
    return ret


def fuse_verts(verts, percentage=0.05):
    r"""
    Work out the span of the points and therefore the range for fusing them together
    then call fuse
    """
    # Work out largest span
    x_span = max(verts[:, 0]) - min(verts[:, 0])
    y_span = max(verts[:, 1]) - min(verts[:, 1])
    if x_span > y_span:
        tolerance = x_span*percentage
    else:
        tolerance = y_span*percentage
    # Fuse vertices lying within 5% of the largest span
    return fuse(verts, tolerance)


def PolyArea2D(pts):
    r"""
    returns the area of a 2D polygon given the set of points defining the convex hull
    in correct order
    """
    lines = np.hstack([pts, np.roll(pts, -1, axis=0)])
    area = 0.5*abs(sum(x1*y2-x2*y1 for x1, y1, x2, y2 in lines))
    return area


def PolyPerimeter2D(pts):
    r"""
    returns the perimeter of a 2D polygon given the set of points defining the convex
    hull in correct order
    """
    lines = np.hstack([pts, np.roll(pts, -1, axis=0)])
    perimeter = sum(np.sqrt((x2-x1)**2+(y2-y1)**2) for x1, y1, x2, y2 in lines)
    return perimeter


def PolyWeightedCentroid2D(pts):
    r"""
    returns the perimeter of a 2D polygon given the set of points defining the convex
    hull in correct order
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


def symmetric_difference(list_a, list_b):
    r"""
    Return the combination of two lists without common elements (necessary as sets
    cannot contain mutable objects)
    """
    sym_diff = []
    sorted_list_a = np.sort(list_a)
    sorted_list_b = np.sort(list_b)
    # Add elements in list a if not in list b
    for element_a in sorted_list_a:
        match = False
        for element_b in sorted_list_b:
            if all(element_a == element_b):
                match = True
        if match is False:
            sym_diff.append(element_a)
    # Add elements in list b if not in list a
    for element_b in sorted_list_b:
        match = False
        for element_a in sorted_list_a:
            if all(element_a == element_b):
                match = True
        if match is False:
            sym_diff.append(element_b)
    return sym_diff


def new_point(pairs):
    r"""
    Passed 2 pairs of points defining lines either side of overlapped offset vertices
    need to calculate the new point to offset from given the orientation of the outer
    fibres
    """
    m1, c1 = line_equation(pairs[0])
    m2, c2 = line_equation(pairs[1])

    if m1 == np.inf:
        # line1 is a straight line x = c1
        x = c1
        y = (m2*c1) + c2
    elif (m2 == np.inf):
        # line2 is a straight line x = c2
        x = c2
        y = (m1*c2) + c1
    else:
        try:
            x = (c2-c1) / (m1-m2)
            y = (m1*c2 - m2*c1) / (m1-m2)
        except RuntimeWarning:
            x = 0
            y = 0
    return x, y


def line_equation(points):
    r"""
    Return the gradient and y intercept of a straight line given 2 points
    """
    x_coords, y_coords = zip(*points)
    dy = y_coords[1]-y_coords[0]
    dx = x_coords[1]-x_coords[0]
    if dx == 0:
        m = np.inf
        c = x_coords[1]
    else:
        m = dy / dx
        c = y_coords[1] - m*x_coords[1]

    return m, c


def set_overlap(verts, offset):
    r"""
    Given a set of vertices and a set of offset vertices, evaluate whether any of the
    offset vertices overlap. This is then used to recalculate points from which to
    offset
    """
    dim = len(verts)
    overlap_array = np.zeros(shape=(dim, dim))
    sweep_radius = np.zeros(len(verts))
    for i, zone_centre in enumerate(verts):
        sweep_radius[i] = np.sqrt(dist2(zone_centre, offset[i]))
        for j, test_point in enumerate(offset):
            test_radius = np.sqrt(dist2(zone_centre, test_point))
            if (test_radius < sweep_radius[i]):
                # Fill in both so that they are both recalculated later i overlapping
                # j doesn't necessarily mean j overlaps i
                overlap_array[i][j] = 1
                overlap_array[j][i] = 1
    # Join up overlapping regions of points
    for i in range(dim):
        for j in range(dim):
            # If an overlap exist look at what others exist for that vertex
            if overlap_array[i][j] == 1:
                for k in range(dim):
                    if overlap_array[j][k] == 1 and k != i:
                        overlap_array[i][k] = 1

    line_points_out = line_points(overlap_array)

    return overlap_array, sweep_radius, line_points_out


def line_points(array):
    r"""
    We are passed a square array containing a list of overlap results. rows represent
    vertices and columns represent offset vertices. If an overlap occurs in the span
    between offset j and offset i then a 1 will result in [i][j]. As the vertices are
    in hull order and our aim is to create a new point from which to offset using
    connected fibres we want to identify the correct points to use to define our
    lines.

      e__________ d
       |        | c
       |       /
       |     /
       |   /
       |_/
       a b
    If we have 5 points in a hull a,b,c,d,e where a overlaps with b and c overlaps
    with d (and visa versa) the array will look like:
    [0,1,0,0,0]
    [1,0,0,0,0]
    [0,0,0,1,0]
    [0,0,1,0,0]
    [0,0,0,0,0]

    This means that c and d are within a fibre's width of each other and the line
    between them does not represent the fibre. Instead we want to extend the outer
    lines (bc and de) to see where they would meet and offset from this point. Roll
    up and down to find the first unoverlapped index from which to start each line
    from then go back one to get the two points to form a line.
    """
    dim = np.shape(array)[0]
    index = range(dim)
    line_points = np.ndarray(shape=[dim, dim], dtype=object)
    for i in range(dim):
        for j in range(dim):
            if array[i][j] == 1 and array[j][i] == 1:
                # Roll forwards to find the first unoverlapped index
                k = 1
                while k < dim:
                    if np.roll(array[i], -k, axis=0)[j] == 0:
                        break
                    else:
                        k += 1
                # Save the indices of the first unoverlapped index and the one before
                # to create the line
                forward_line = [np.roll(index, -k, axis=0)[j],
                                np.roll(index, -(k-1), axis=0)[j]]
                forward_line.sort()
                # Roll backwards to find the first unoverlapped index
                k = 1
                while k < dim:
                    if np.roll(array[i], k, axis=0)[j] == 0:
                        break
                    else:
                        k += 1
                # Save the indices of the first unoverlapped index and the one before
                # to create the line
                backward_line = [np.roll(index, k, axis=0)[j],
                                 np.roll(index, (k - 1), axis=0)[j]]
                backward_line.sort()
                line_points[i][j] = (forward_line, backward_line)

    return line_points


def all_overlap(array):
    r"""
    Find out whether all offset vertices (columns) are overlapped by at least one
    other. If so then throat is fully occluded
    """
    dim = np.shape(array)[0]
    overlap = [False]*dim
    all_overlap = False
    for i in range(dim):
        for j in range(dim):
            if array[j][i] == 1:
                overlap[i] = True
    if sum(overlap) == dim:
        all_overlap = True

    return all_overlap


def scale(network, scale_factor=[1, 1, 1], preserve_vol=True):
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
    >>> vo.scale(network=pn,scale_factor=[2,1,1])
    >>> Vol2 = vo.vertex_dimension(pn,B1,B2)
    >>> np.around(Vol-Vol2,5) == 0.0
    True
    >>> vo.scale(network=pn,scale_factor=[2,1,1],preserve_vol=False)
    >>> Vol3 = vo.vertex_dimension(pn,B1,B2)
    >>> np.around(Vol3/Vol,5)
    2.0

    """
    from scipy.special import cbrt
    import scipy as sp
    scale_factor = np.asarray(scale_factor)
    if preserve_vol is True:
        scale_factor = scale_factor/(cbrt(sp.prod(scale_factor)))
    network['pore.coords'] = network['pore.coords']*scale_factor
    # Cycle through all vertices of all pores updating vertex values
    for pore in network.pores():
        for i, vert in network['pore.vert_index'][pore].items():
            network['pore.vert_index'][pore][i] = \
                network['pore.vert_index'][pore][i]*scale_factor
    # Cycle through all vertices of all throats updating vertex values
    for throat in network.throats():
        for i, vert in network['throat.vert_index'][throat].items():
            network['throat.vert_index'][throat][i] = \
                network['throat.vert_index'][throat][i]*scale_factor


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
        pore_vol = sum(network['pore.volume'])
    except KeyError:
        print('Geometries must be assigned first')
        pore_vol = 0
    porosity = pore_vol/domain_vol
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
        print('Something is wrong: ' + str(tot_angle))

    return 1 / np.cos(np.array([theta_x, theta_y, theta_z]))


def print_throat(geom, throats_in):
    r"""
    Print a given throat or list of throats accepted as [1, 2, 3, ..., n]
    e.g geom.print_throat([34, 65, 99])
    Original vertices plus offset vertices are rotated to align with
    the z-axis and then printed in 2D
    """
    import matplotlib.pyplot as plt
    throats = []
    for throat in throats_in:
        if throat in range(geom.num_throats()):
            throats.append(throat)
        else:
            print('Throat: ' + str(throat) + ' not part of geometry')
    if len(throats) > 0:
        verts = geom['throat.vertices'][throats]
        offsets = geom['throat.offset_vertices'][throats]
        normals = geom['throat.normal'][throats]
        coms = geom['throat.centroid'][throats]
        incentre = geom['throat.incentre'][throats]
        inradius = 0.5*geom['throat.indiameter'][throats]
        for i in range(len(throats)):
            fig = plt.figure()
            vert_2D = tr.rotate_and_chop(verts[i], normals[i], [0, 0, 1])
            hull = ConvexHull(vert_2D, qhull_options='QJ Pp')
            for simplex in hull.simplices:
                plt.plot(vert_2D[simplex, 0], vert_2D[simplex, 1], 'k-', linewidth=2)
            plt.scatter(vert_2D[:, 0], vert_2D[:, 1])
            offset_2D = tr.rotate_and_chop(offsets[i], normals[i], [0, 0, 1])
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
            centroid = tr.rotate_and_chop(coms[i], normals[i], [0, 0, 1])
            incent = tr.rotate_and_chop(incentre[i], normals[i], [0, 0, 1])
            plt.scatter(centroid[0][0], centroid[0][1])
            # Plot incircle
            t = np.linspace(0, 2*np.pi, 200)
            u = inradius[i]*np.cos(t)+incent[0][0]
            v = inradius[i]*np.sin(t)+incent[0][1]
            plt.plot(u, v, 'r-')
            fig.show()
    else:
        print('Please provide throat indices')


def print_pore(geom, pores, fig=None, axis_bounds=None):
    r"""
    Print all throats around a given pore or list of pores accepted
    as [1, 2, 3, ..., n]
    e.g geom.print_pore([34, 65, 99])
    Original vertices plus offset vertices used to create faces and
    then printed in 3D
    To print all pores (n)
    pore_range = np.arange(0,n-1,1)
    geom.print_pore(pore_range)
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    return_fig = False
    if len(pores) > 0:
        net_pores = geom.map_pores(geom._net, pores)
        centroids = geom['pore.centroid'][pores]
        throats = geom._net.map_throats(geom,
                                        net_throats,
                                        return_mapping=True)['target']
        # Can't create volume from one throat
        if 1 <= len(throats):
            verts = geom['throat.vertices'][throats]
            normals = geom['throat.normal'][throats]
            # Get verts in hull order
            ordered_verts = []
            for i in range(len(verts)):
                vert_2D = tr.rotate_and_chop(verts[i], normals[i], [0, 0, 1])
                hull = ConvexHull(vert_2D, qhull_options='QJ Pp')
                ordered_verts.append(verts[i][hull.vertices])
            offsets = geom['throat.offset_vertices'][throats]
            ordered_offs = []
            for i in range(len(offsets)):
                offs_2D = tr.rotate_and_chop(offsets[i], normals[i], [0, 0, 1])
                offs_hull = ConvexHull(offs_2D, qhull_options='QJ Pp')
                ordered_offs.append(offsets[i][offs_hull.vertices])
            # Get domain extents for setting axis
            if axis_bounds is None:
                [xmin, xmax, ymin, ymax, zmin, zmax] = \
                    vertex_dimension(geom._net, pores, parm='minmax')
            else:
                [xmin, xmax, ymin, ymax, zmin, zmax] = axis_bounds
            if fig is None:
                fig = plt.figure()
            else:
                return_fig is True
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
            ax.scatter(centroids[:, 0], centroids[:, 1], centroids[:, 2], c='y')
            plt.show()
        else:
            print_throat(throats)
    else:
        print('Please provide pore indices')
    if return_fig is True:
        return fig
