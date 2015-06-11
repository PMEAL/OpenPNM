r"""
===============================================================================
pore_centroid -- N.B Neither method is accurate but is quicker than that
                 implemented in pore.volume
===============================================================================

"""
import scipy as _sp


def voronoi(geometry, pore_vertices='pore.vertices', **kwargs):
    r"""
    Calculate the centroid of the pore from the voronoi vertices - C.O.M
    """
    verts = geometry[pore_vertices]
    value = _sp.ndarray([len(verts), 3])
    for i, vert in enumerate(verts):
        value[i] = _sp.array([vert[:, 0].mean(),
                              vert[:, 1].mean(),
                              vert[:, 2].mean()])
    return value


def voronoi2(geometry, vertices='throat.centroid', **kwargs):
    r"""
    Calculate the centroid from the mean of the throat centroids
    """
    value = _sp.ndarray([geometry.num_pores(), 3])
    pore_map = geometry.map_pores(geometry.pores(), geometry._net)
    for geom_pore, net_pore in pore_map:
        net_throats = geometry._net.find_neighbor_throats(net_pore)
        geom_throats = geometry._net.map_throats(net_throats, geometry)[:, 1]
        verts = geometry[vertices][geom_throats]
        value[geom_pore] = _sp.array([verts[:, 0].mean(),
                                      verts[:, 1].mean(),
                                      verts[:, 2].mean()])
    return value
