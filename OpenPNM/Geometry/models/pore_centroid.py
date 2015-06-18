r"""
===============================================================================
pore_centroid
===============================================================================

"""
import scipy as _sp


def voronoi(network, geometry, vertices='throat.centroid', **kwargs):
    r"""
    Calculate the centroid from the mean of the throat centroids
    """
    value = _sp.ndarray([geometry.num_pores(), 3])
    pore_map = geometry.map_pores(target=network,
                                  pores=geometry.pores(),
                                  return_mapping=True)
    for i, net_pore in enumerate(pore_map['target']):
        geom_pore = pore_map['source'][i]
        net_throats = network.find_neighbor_throats(net_pore)
        geom_throats = network.map_throats(target=geometry,
                                           throats=net_throats,
                                           return_mapping=False)
        verts = geometry[vertices][geom_throats]
        " Ignore all zero centroids "
        verts = verts[~_sp.all(verts == 0, axis=1)]
        value[geom_pore] = _sp.mean(verts, axis=0)

    return value
