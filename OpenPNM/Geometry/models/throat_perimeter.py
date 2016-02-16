r"""
===============================================================================
Submodule -- throat_perimeter
===============================================================================

"""
import scipy as _sp
import OpenPNM.Utilities.vertexops as vo
from scipy.spatial import ConvexHull


def cylinder(geometry, throat_diameter='throat.diameter', **kwargs):
    r"""
    Calcuate the throat perimeter assuming a circular cross-section

    Parameters
    ----------
    geometry : OpenPNM Geometry object
        The Geometry object with which this model is associated.  This is
        needed to access the ``throat_diameter`` values.

    throat_diameter : string
        The dictionary key of the array containing the throat diameter values
    """
    return geometry[throat_diameter]*_sp.constants.pi


def cuboid(geometry, throat_diameter='throat.diameter', **kwargs):
    r"""
    Calcuate the throat perimeter assuming a square cross-section

    Parameters
    ----------
    geometry : OpenPNM Geometry object
        The Geometry object with which this model is associated.  This is
        needed to access the ``throat_diameter`` values.

    throat_diameter : string
        The dictionary key of the array containing the throat diameter values
    """
    return geometry[throat_diameter]*4


def voronoi(geometry, **kwargs):
    r"""
    Uses the Voronoi vertices and throat normals to work out the perimeter
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
