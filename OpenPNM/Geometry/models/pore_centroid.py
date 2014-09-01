r"""
===============================================================================
Submodule -- pore_centroid
===============================================================================

"""
import scipy as _sp

def voronoi(geometry,
            pore_vertices='pore.vertices',
            **kwargs):
    r"""
    Calculate the centroid of the pore from the voronoi vertices - C.O.M
    """
    network = geometry._net    
    pores = network.pores()
    verts = network[pore_vertices][pores]
    value = _sp.ndarray(len(verts),dtype=object)
    for i,vert in enumerate(verts):
        value[i] = _sp.array([vert[:,0].mean(),vert[:,1].mean(),vert[:,2].mean()])
    return value