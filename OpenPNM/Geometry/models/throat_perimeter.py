r"""
===============================================================================
Submodule -- throat_perimeter
===============================================================================

"""
import scipy as _sp

def voronoi(geometry,
            **kwargs):
    r"""
    Use the Voronoi verts and throat normals to work out the perimeter
    """
    Nt = geometry.num_throats()    
    verts = geometry['throat.offset_vertices']
    normals = geometry['throat.normal']
    perimeter = _sp.ndarray(Nt)
    for i in range(Nt):
        verts_2D = geometry._rotate_and_chop(verts[i],normals[i],[0,0,1])
        perimeter[i] = geometry._PolyPerimeter2D(verts_2D)
    
    return perimeter