r"""
===============================================================================
throat_vertices -- update throat vertices from Vornoi object
===============================================================================

"""
import scipy as _sp

def voronoi(network,
            geometry,
              **kwargs):
    r"""
    Update the pore vertices from the voronoi vertices
    """    
    throats = geometry["throat.map"]    
    value = _sp.ndarray(len(throats),dtype=object)
    for i in range(len(throats)):
        value[i]=network._vor.vertices[network["throat.vert_index"][throats[i]]]
    return value