r"""
===============================================================================
Submodule -- throat_area
===============================================================================

"""
import scipy as _sp

def cylinder(geometry,
             throat_diameter='throat.diameter',
             **kwargs):
    r"""
    Calculate throat area for a cylindrical throat
    """
    diams = geometry[throat_diameter]
    value = _sp.constants.pi/4*(diams)**2
    return value

def cuboid(geometry,
           throat_diameter='throat.diameter',
           **kwargs):
    r"""
    Calculate throat area for a cuboid throat
    """
    diams = geometry[throat_diameter]
    value = (diams)**2
    return value
        
def voronoi(geometry,
            **kwargs):
    r"""
    Use the Voronoi verts and throat normals to work out the area
    """
    Nt = geometry.num_throats()    
    verts = geometry['throat.offset_vertices']
    normals = geometry['throat.normal']
    area = _sp.ndarray(Nt)
    for i in range(Nt):
        verts_2D = geometry._rotate_and_chop(verts[i],normals[i],[0,0,1])
        area[i] = geometry._PolyArea2D(verts_2D)
    
    return area