r"""
===============================================================================
throat_vertices -- update throat vertices from Vornoi object
===============================================================================

"""
import scipy as _sp


def voronoi(network, geometry, **kwargs):
    r"""
    Update the pore vertices from the voronoi vertices
    """
    throats = geometry.map_throats(network, geometry. throats())
    value = _sp.ndarray(len(throats), dtype=object)
    for i in range(len(throats)):
        value[i] = \
            _sp.asarray(list(network['throat.vert_index'][throats[i]].values()))
    return value
