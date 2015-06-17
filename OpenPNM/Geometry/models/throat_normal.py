r"""
===============================================================================
throat_normal -- calculate from the normal vector to the throat vertices
===============================================================================

"""
import scipy as sp
from scipy.spatial import ConvexHull


def voronoi(network, geometry, **kwargs):
    r"""
    Update the throat normals from the voronoi vertices
    """
    verts = geometry["throat.vertices"]
    value = sp.ndarray([len(verts), 3])
    for i in range(len(verts)):
        if len(sp.unique(verts[i][:, 0])) == 1:
            verts_2d = sp.vstack((verts[i][:, 1], verts[i][:, 2])).T
        elif len(sp.unique(verts[i][:, 1])) == 1:
            verts_2d = sp.vstack((verts[i][:, 0], verts[i][:, 2])).T
        else:
            verts_2d = sp.vstack((verts[i][:, 0], verts[i][:, 1])).T
        hull = ConvexHull(verts_2d, qhull_options='QJ Pp')
        sorted_verts = verts[i][hull.vertices]
        v1 = sorted_verts[1]-sorted_verts[0]
        v2 = sorted_verts[-1]-sorted_verts[0]
        value[i] = sp.cross(v1, v2)

    return value
