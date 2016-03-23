r"""
===============================================================================
Submodule -- throat_centroid
===============================================================================

"""
import scipy as _sp
from transforms3d import _gohlketransforms as tr
from OpenPNM.Utilities import vertexops as vo
from scipy.spatial import ConvexHull


def centre_of_mass(geometry, vertices='throat.offset_vertices', **kwargs):
    r"""
    Calculate the centre of mass of the throat from the voronoi vertices.
    """
    Nt = geometry.num_throats()
    outer_verts = geometry['throat.vertices']
    offset_verts = geometry[vertices]
    normal = geometry['throat.normal']
    z_axis = [0, 0, 1]
    value = _sp.ndarray([Nt, 3])
    for i in range(Nt):
        if len(offset_verts[i]) > 2:
            verts = offset_verts[i]
        elif len(outer_verts[i]) > 2:
            verts = outer_verts[i]
        else:
            verts = []
        if len(verts) > 0:
            # For boundaries some facets will already be aligned with the axis -
            # if this is the case a rotation is unnecessary and could also cause
            # problems
            angle = tr.angle_between_vectors(normal[i], z_axis)
            if angle == 0.0 or angle == _sp.pi:
                "We are already aligned"
                rotate_input = False
                facet = verts
            else:
                rotate_input = True
                M = tr.rotation_matrix(tr.angle_between_vectors(normal[i], z_axis),
                                       tr.vector_product(normal[i], z_axis))
                facet = _sp.dot(verts, M[:3, :3].T)
            # Now we have a rotated facet aligned with the z axis - make 2D
            facet_2D = _sp.column_stack((facet[:, 0], facet[:, 1]))
            z = _sp.unique(_sp.around(facet[:, 2], 10))
            if len(z) == 1:
                # We need the vertices arranged in order so perform a convex hull
                hull = ConvexHull(facet_2D)
                ordered_facet_2D = facet_2D[hull.vertices]
                # Call the routine to calculate an area wighted centroid from the
                # 2D polygon
                COM_2D = vo.PolyWeightedCentroid2D(ordered_facet_2D)
                COM_3D = _sp.hstack((COM_2D, z))
                # If we performed a rotation we need to rotate back
                if (rotate_input):
                    MI = tr.inverse_matrix(M)
                    # Unrotate the offset coordinates using the inverse of the
                    # original rotation matrix
                    value[i] = _sp.dot(COM_3D, MI[:3, :3].T)
                else:
                    value[i] = COM_3D
            else:
                logger.error('Rotation Failed: ' + str(_sp.unique(facet[:, 2])))

    return value
