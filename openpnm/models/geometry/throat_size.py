r"""
"""
from openpnm.models import misc as _misc
import numpy as _np


def weibull(target, shape, scale, loc, seeds='throat.seed'):
    return _misc.weibull(target=target, shape=shape, scale=scale, loc=loc,
                         seeds=seeds)


weibull.__doc__ = _misc.weibull.__doc__


def normal(target, scale, loc, seeds='throat.seed'):
    return _misc.normal(target=target, scale=scale, loc=loc, seeds=seeds)


normal.__doc__ = _misc.normal.__doc__


def generic_distribution(target, func, seeds='throat.seed'):
    return _misc.generic_distribution(target=target, func=func, seeds=seeds)


generic_distribution.__doc__ = _misc.generic_distribution.__doc__


def random(target, seed=None, num_range=[0, 1]):
    return _misc.random(target=target, element='throat', seed=seed,
                        num_range=num_range)


random.__doc__ = _misc.random.__doc__


def from_neighbor_pores(target, prop='pore.diameter', mode='min'):
    return _misc.from_neighbor_pores(target=target, prop=prop,
                                     mode=mode)


from_neighbor_pores.__doc__ = _misc.from_neighbor_pores.__doc__


def equivalent_diameter(target, throat_area='throat.area',
                        throat_shape='circle'):
    r"""
    Calculates the diameter of a cirlce or edge-length of a sqaure with same
    area as the throat.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    thorat_area : string
        The dictionary key to the throat area values

    throat_shape : string
        The shape cross-sectional shape of the throat to assume when
        back-calculating from the area.  Options are 'circle' (default) or
        'square'.

    Returns
    -------
    value : NumPy ndarray
        Array containing throat equivalent diameter.

    """
    area = target[throat_area]
    if throat_shape.startswith('circ'):
        value = 2*_np.sqrt(area/_np.pi)
    elif throat_shape.startswith('square'):
        value = _np.sqrt(area)
    return value
