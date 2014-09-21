r"""
===============================================================================
throat_offset_vertices -- Offeset throat vertices using a fibre radius parameter
===============================================================================

"""
import scipy as sp

def voronoi(network,
            geometry,
            offset,
              **kwargs):
    r"""
    Offset the throat vertices effectively erroding the facet by the offset distance supplied
    """    
    Nt = geometry.num_throats()
    area = sp.ndarray(Nt)
    perimeter = sp.ndarray(Nt)
    offset_verts = sp.ndarray(Nt,dtype=object)
    offset_error = sp.ndarray(Nt)
    for i in range(Nt):            
        throat_verts=geometry["throat.vertices"][i]
        throat_normal=geometry["throat.normal"][i]
        area[i],perimeter[i],offset_verts[i],offset_error[i] = geometry._get_throat_geom(throat_verts,throat_normal,offset)    
    
    return offset_verts