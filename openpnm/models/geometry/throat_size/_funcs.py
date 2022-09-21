import numpy as _np
from openpnm.models import misc as _misc
from openpnm.models.geometry import _geodocs


__all__ = ["weibull",
           "normal",
           "generic_distribution",
           "random",
           "from_neighbor_pores",
           "equivalent_diameter"]


def weibull(network, shape, scale, loc, seeds='throat.seed'):
    return _misc.weibull(network, shape=shape, scale=scale, loc=loc,
                         seeds=seeds)


weibull.__doc__ = _misc.weibull.__doc__


def normal(network, scale, loc, seeds='throat.seed'):
    return _misc.normal(network, scale=scale, loc=loc, seeds=seeds)


normal.__doc__ = _misc.normal.__doc__


def generic_distribution(network, func, seeds='throat.seed'):
    return _misc.generic_distribution(network, func=func, seeds=seeds)


generic_distribution.__doc__ = _misc.generic_distribution.__doc__


def random(network, seed=None, num_range=[0, 1]):
    return _misc.random(network, element='throat', seed=seed,
                        num_range=num_range)


random.__doc__ = _misc.random.__doc__


def from_neighbor_pores(network, prop='pore.diameter', mode='min'):
    return _misc.from_neighbor_pores(network, prop=prop,
                                     mode=mode)


from_neighbor_pores.__doc__ = _misc.from_neighbor_pores.__doc__


@_geodocs
def equivalent_diameter(
    network,
    throat_area='throat.cross_sectional_area',
    throat_shape='circle',
):
    r"""
    Calculates the diameter of a cirlce or edge-length of a sqaure with same
    area as the throat.

    Parameters
    ----------
    %(network)s
    %(At)s
    throat_shape : str
        The shape cross-sectional shape of the throat to assume when
        back-calculating from the area.  Options are 'circle' (default) or
        'square'.

    Returns
    -------
    diameters : ndarray
        A numpy ndarray containing throat diameter values

    """
    area = network[throat_area]
    if throat_shape.startswith('circ'):
        value = 2*_np.sqrt(area/_np.pi)
    elif throat_shape.startswith('square'):
        value = _np.sqrt(area)
    return value
