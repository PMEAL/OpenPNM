r"""
===============================================================================
Submodule -- throat_centroid
===============================================================================

"""
import scipy as _sp

def voronoi(geometry,
            **kwargs):
    r"""
    Calculate the centroid of the throat from the voronoi vertices - C.O.M
    """
    network = geometry._net    
    throats = network.throats(geometry.name)
    verts = network['throat.offset_verts'][throats]
    value = _sp.ndarray(len(verts),dtype=object)
    for i,vert in enumerate(verts):
        value[i] = _sp.array([vert[:,0].mean(),vert[:,1].mean(),vert[:,2].mean()])
    return value