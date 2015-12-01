r"""
===============================================================================
Submodule -- throat_perimeter
===============================================================================

"""
import scipy as _sp
import OpenPNM.Utilities.vertexops as vo
from scipy.spatial import ConvexHull


def voronoi(geometry, **kwargs):
    r"""
    Use the Voronoi verts and throat normals to work out the perimeter
    """
    Nt = geometry.num_throats()
    verts = geometry['throat.offset_vertices']
    normals = geometry['throat.normal']
    perimeter = _sp.ndarray(Nt)
    for i in range(Nt):
        if len(verts[i]) > 2:
            verts_2D = vo.rotate_and_chop(verts[i], normals[i], [0, 0, 1])
            # Get in hull order
            hull = ConvexHull(verts_2D, qhull_options='QJ Pp')
            verts_2D = verts_2D[hull.vertices]
            perimeter[i] = vo.PolyPerimeter2D(verts_2D)
        else:
            perimeter[i] = 0.0
    return perimeter
