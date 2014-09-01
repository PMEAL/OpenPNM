r"""
===============================================================================
Submodule -- throat_centroid
===============================================================================

"""
import scipy as _sp

def voronoi(network,
            **kwargs):
    r"""
    Calculate the centroid of the pore from the voronoi vertices - C.O.M
    """
    throats = network.throats()
    verts = network['throat.offset_verts'][throats]
    value = _sp.ndarray(len(verts),dtype=object)
    for i,vert in enumerate(verts):
        value[i] = _sp.array([vert[:,0].mean(),vert[:,1].mean(),vert[:,2].mean()])
    return value