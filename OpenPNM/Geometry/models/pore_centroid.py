r"""
===============================================================================
pore_centroid
===============================================================================

"""
import scipy as _sp


def voronoi(geometry, vertices='throat.centroid', **kwargs):
    r"""
    Calculate the centroid from the mean of the throat centroids

    Parameters
    ----------
    geometry : OpenPNM Geometry object
        The Geometry object which this model is associated.  This is needed to
        access the values of the ``vertices``.

    vertices : string
        The dictionary contain of the array containing the throat vertice
        coordiantes.  The default is 'throat.centroid'

    """
    network = geometry._net
    value = _sp.ndarray([geometry.num_pores(), 3])
    value.fill(0.0)
    pore_map = geometry.map_pores(target=network,
                                  pores=geometry.pores(),
                                  return_mapping=True)
    for i, net_pore in enumerate(pore_map['target']):
        geom_pore = pore_map['source'][i]
        net_throats = network.find_neighbor_throats(net_pore)
        geom_throats = network.map_throats(target=geometry,
                                           throats=net_throats,
                                           return_mapping=True)['target']
        verts = geometry[vertices][geom_throats]
        " Ignore all zero centroids "
        verts = verts[~_sp.all(verts == 0, axis=1)]
        if len(verts) > 0:
            value[geom_pore] = _sp.mean(verts, axis=0)

    return value
