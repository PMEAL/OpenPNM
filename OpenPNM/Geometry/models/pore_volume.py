r"""
===============================================================================
pore_volume --
===============================================================================

"""
import scipy as _sp
import numpy as np
from scipy.spatial import Delaunay
import OpenPNM.Utilities.misc as misc


def sphere(geometry, pore_diameter='pore.diameter', **kwargs):
    r"""
    Calculate pore volume from diameter assuming a spherical pore body

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key to the pore diameter values

    """
    diams = geometry[pore_diameter]
    value = _sp.pi/6*diams**3
    return value


def cube(geometry, pore_diameter='pore.diameter', **kwargs):
    r"""
    Calculate pore volume from diameter assuming a cubic pore body

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key to the pore diameter values

    """
    diams = geometry[pore_diameter]
    value = diams**3
    return value


def _get_hull_volume(points):
    r"""
    Calculate the volume of a set of points by dividing the bounding surface
    into triangles and working out the volume of all the pyramid elements connected
    to the volume centroid
    """
    # Remove any duplicate points - this messes up the triangulation
    points = _sp.asarray(misc.unique_list(np.around(points, 10)))
    try:
        tri = Delaunay(points, qhull_options='QJ Pp')
    except _sp.spatial.qhull.QhullError:
        print(points)
    # We only want points included in the convex hull to calculate the centroid
    hull_centroid = _sp.array([points[:, 0].mean(),
                               points[:, 1].mean(),
                               points[:, 2].mean()])
    hull_volume = 0.0
    pyramid_COMs = []
    for ia, ib, ic in tri.convex_hull:
        # Points making each triangular face
        # Collection of co-ordinates of each point in this face
        face_x = points[[ia, ib, ic]][:, 0]
        face_y = points[[ia, ib, ic]][:, 1]
        face_z = points[[ia, ib, ic]][:, 2]
        # Average of each co-ordinate is the centroid of the face
        face_centroid = [face_x.mean(), face_y.mean(), face_z.mean()]
        face_centroid_vector = face_centroid - hull_centroid
        # Vectors of the sides of the face used to find normal vector and area
        vab = points[ib] - points[ia]
        vac = points[ic] - points[ia]
        vbc = points[ic] - points[ib]
        # As vectors are co-planar the cross-product will produce the normal
        # vector of the face
        face_normal = _sp.cross(vab, vac)
        try:
            face_unit_normal = face_normal/_sp.linalg.norm(face_normal)
        except RuntimeWarning:
            print('Pore Volume Error:' + str(vab) + ' ' + str(vac))
        # As triangles are orientated randomly in 3D we could either transform
        # co-ordinates to align with a plane and perform 2D operations to work out
        # the area or we could work out the lengths of each side and use Heron's
        # formula which is easier. Using Delaunay traingulation will always produce
        # triangular faces but if dealing with other polygons co-ordinate transfer
        # may be necessary
        a = _sp.linalg.norm(vab)
        b = _sp.linalg.norm(vbc)
        c = _sp.linalg.norm(vac)
        # Semiperimeter
        s = 0.5*(a + b + c)
        face_area = _sp.sqrt(s*(s-a)*(s-b)*(s-c))
        # Now the volume of the pyramid section defined by the 3 face points and the
        # hull centroid can be calculated
        pyramid_volume = _sp.absolute(_sp.dot(face_centroid_vector,
                                              face_unit_normal)*face_area/3)
        # Each pyramid is summed together to calculate the total volume
        hull_volume += pyramid_volume
        # The Centre of Mass will not be the same as the geometrical centroid. A
        # weighted adjustment can be calculated from the pyramid centroid and volume
        vha = points[ia] - hull_centroid
        vhb = points[ib] - hull_centroid
        vhc = points[ic] - hull_centroid
        pCOM = ((vha+vhb+vhc)/4)*pyramid_volume
        pyramid_COMs.append(pCOM)
    if _sp.isnan(hull_volume):
        hull_volume = 0.0
    if hull_volume > 0:
        hull_COM = hull_centroid + _sp.mean(_sp.asarray(pyramid_COMs),
                                            axis=0) / hull_volume
    else:
        hull_COM = hull_centroid

    return hull_volume, hull_COM


def voronoi(network, geometry, **kwargs):
    r"""
    Calculate volume from the convex hull of the offset vertices making the throats
    surrounding the pore. Also calculate the centre of mass for the volume
    """
    pores = geometry.map_pores(network, geometry.pores())
    Np = len(pores)
    volume = _sp.zeros(Np)
    com = _sp.zeros([Np, 3])
    for i in range(Np):
        throat_vert_list = []
        net_throats = network.find_neighbor_throats([pores[i]])
        geom_throats = network.map_throats(target=geometry,
                                           throats=net_throats,
                                           return_mapping=True)['target']
        if len(geom_throats) > 1:
            for throat in geom_throats:
                geom_throat_verts = geometry["throat.offset_vertices"][throat]
                if geom_throat_verts is not None:
                    for j in range(len(geom_throat_verts)):
                        throat_vert_list.append(geom_throat_verts[j])
            throat_array = _sp.asarray(throat_vert_list)
            if len(throat_array) > 4:
                volume[i], com[i] = _get_hull_volume(throat_array)
            else:
                volume[i] = 0
        elif len(geom_throats) == 1 and 'throat.centroid' in geometry.props():
                com[i] = geometry['throat.centroid'][geom_throats]
                volume[i] = 0
    # Find any pores with centroids at origin and use the mean of the pore vertices
    # instead. Not doing this messes up hydraulic conductances using centre to centre
    ps = np.where(~com.any(axis=1))[0]
    if len(ps) > 0:
        for pore in ps:
            com[pore] = np.mean(geometry['pore.vertices'][pore], axis=0)
    geometry['pore.centroid'] = com

    return volume
