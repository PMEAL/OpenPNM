r"""
===============================================================================
Submodule -- throat_area
===============================================================================

"""
import scipy as _sp
import OpenPNM.Utilities.transformations as tr
import OpenPNM.Utilities.vertexops as vo


def cylinder(geometry, throat_diameter='throat.diameter', **kwargs):
    r"""
    Calculate throat area for a cylindrical throat
    """
    diams = geometry[throat_diameter]
    value = _sp.constants.pi/4*(diams)**2
    return value


def cuboid(geometry, throat_diameter='throat.diameter', **kwargs):
    r"""
    Calculate throat area for a cuboid throat
    """
    diams = geometry[throat_diameter]
    value = (diams)**2
    return value


def voronoi(geometry, **kwargs):
    r"""
    Use the Voronoi verts and throat normals to work out the area
    """
    Nt = geometry.num_throats()
    verts = geometry['throat.offset_vertices']
    normals = geometry['throat.normal']
    area = _sp.ndarray(Nt)
    for i in range(Nt):
        if len(verts[i]) > 2:
            verts_2D = tr.rotate_and_chop(verts[i], normals[i], [0, 0, 1])
            area[i] = vo.PolyArea2D(verts_2D)
        else:
            area[i] = 0.0

    return area
