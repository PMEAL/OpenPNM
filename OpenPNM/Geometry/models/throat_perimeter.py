r"""
===============================================================================
Submodule -- throat_perimeter
===============================================================================

"""
import scipy as _sp

def voronoi(geometry,
            fibre_rad,
            throat_area='throat.perimeter',
            **kwargs):
    r"""
    Use the Voronoi verts and throat normals to work out the perimeter
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
    return perimeter