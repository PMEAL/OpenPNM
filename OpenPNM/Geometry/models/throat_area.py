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
            fibre_rad,
            throat_area='throat.area',
            **kwargs):
    r"""
    Use the Voronoi verts and throat normals to work out the area
    """
    Nt = geometry.num_throats()    
    throats = geometry['throat.map']
    network = geometry._net
    verts = network['throat.verts'][throats]
    normals = network['throat.normals'][throats]
    area = _sp.ndarray(Nt)
    perimeter = _sp.ndarrat(Nt)
    offset_verts = _sp.ndarray(Nt,dtype=object)
    error = _sp.ndarray(Nt)
    for i in _sp.arange(0,Nt):
        area[i],perimeter[i],offset_verts[i],error[i]=geometry._get_throat_geom(verts[i],normals[i],fibre_rad)
    return area