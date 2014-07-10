r"""
===============================================================================
Submodule -- pore_centroid
===============================================================================

"""
import scipy as _sp

def voronoi(network,pores,**kwargs):
    r"""
    Calculate the centroid of the pore from the voronoi vertices - C.O.M
    """
    verts = network['pore.vertices'][pores]
    value = _sp.ndarray(len(verts),dtype=object)
    for i,vert in enumerate(verts):
        value[i] = _sp.array([vert[:,0].mean(),vert[:,1].mean(),vert[:,2].mean()])
    return value