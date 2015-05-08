r"""
===============================================================================
pore_vertices -- update pore vertices from Vornoi object
===============================================================================

"""
import scipy as _sp


def voronoi(network, geometry, **kwargs):
    r"""
    Update the pore vertices from the voronoi vertices
    """
    pores = geometry.map_pores(network, geometry.pores())
    value = _sp.ndarray(len(pores), dtype=object)
    for i in range(len(pores)):
        value[i] = _sp.asarray(list(network['pore.vert_index'][pores[i]].values()))
    return value
