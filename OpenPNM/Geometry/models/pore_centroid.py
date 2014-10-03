r"""
===============================================================================
pore_centroid -- 
===============================================================================

"""
import scipy as _sp

def voronoi(geometry,
            pore_vertices='pore.vertices',
            **kwargs):
    r"""
    Calculate the centroid of the pore from the voronoi vertices - C.O.M
    """
    #network = geometry._net    
    #pores = geometry['pore.map']
    verts = geometry[pore_vertices]
    value = _sp.ndarray([len(verts),3])
    for i,vert in enumerate(verts):
        value[i] = _sp.array([vert[:,0].mean(),vert[:,1].mean(),vert[:,2].mean()])
    return value

def voronoi2(geometry,
             vertices='throat.centroid',
             **kwargs):
    r"""
    Calculate the centroid from the mean of the throat centroids
    """
    value = _sp.ndarray([geometry.num_pores(),3])
    for geom_pore,net_pore in enumerate(geometry["pore.map"]):
        throats = geometry._net.find_neighbor_throats(net_pore)
        verts = geometry[vertices][throats]
        value[geom_pore]=_sp.array([verts[:,0].mean(),verts[:,1].mean(),verts[:,2].mean()])
    return value