from openpnm.models import misc
from .throat_misc import neighbor as _neighbor
import nump as np


def weibull(target, shape, scale, loc, seeds='throat.seed'):
    return misc.weibull(target=target, shape=shape, scale=scale, loc=loc,
                        seeds=seeds)


weibull.__doc__ = misc.weibull.__doc__


def normal(target, scale, loc, seeds='throat.seed'):
    return misc.normal(target=target, scale=scale, loc=loc, seeds=seeds)


normal.__doc__ = misc.normal.__doc__


def generic(target, func, seeds='throat.seed'):
    return misc.generic(target=target, func=func, seeds=seeds)


generic.__doc__ = misc.generic.__doc__


def neighbor(target, pore_prop, mode='min'):
    return _neighbor(target=target, pore_prop=pore_prop, mode=mode)


neighbor.__doc__ = _neighbor.__doc__


def random(target, seed=None, num_range=[0, 1]):
    return misc.random(target=target, element='throat', seed=seed,
                       num_range=num_range)


random.__doc__ = misc.random.__doc__


def equivalent_circle(target, throat_area='throat.area'):
    r"""
    Calculates the diameter of a cirlce with same area as the throat.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    thorat_area : string
        The dictionary key to the throat area values
    """
    areas = target[throat_area]
    value = 2*np.sqrt(areas/np.pi)
    return value
