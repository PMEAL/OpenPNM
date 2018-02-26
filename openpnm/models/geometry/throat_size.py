r"""
===============================================================================
Submodule -- throat_diameter
===============================================================================

"""
from . import misc as _misc
from .throat_misc import neighbor as _neighbor
import scipy as _sp


def weibull(target, shape, scale, loc, seeds='throat.seed'):
    return _misc.weibull(target=target, shape=shape, scale=scale, loc=loc,
                         seeds=seeds)


weibull.__doc__ = _misc.weibull.__doc__


def normal(target, scale, loc, seeds='throat.seed'):
    return _misc.normal(target=target, scale=scale, loc=loc, seeds=seeds)


normal.__doc__ = _misc.normal.__doc__


def generic(target, func, seeds='throat.seed'):
    return _misc.generic(target=target, func=func, seeds=seeds)


generic.__doc__ = _misc.generic.__doc__


def neighbor(target, pore_prop, mode='min'):
    return _neighbor(target=target, pore_prop=pore_prop, mode=mode)


neighbor.__doc__ = _neighbor.__doc__


def random(target, seed=None, num_range=[0, 1]):
    r"""
    Assign throat sizes from a random

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    seed : int
        The starting seed value to send to Scipy's random number generator.
        The default is None, which means different distribution is returned
        each time the model is run.

    num_range : list
        A two element list indicating the low and high end of the returned
        numbers.  The default is [0, 1] but this can be adjusted to produce
        throat sizes directly; for instance pores between 10 and 100 um can be
        generated with ``num_range = [0.00001, 0.0001]``.
    """
    return _misc.random(target=target, element='throat', seed=seed,
                        num_range=num_range)


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
    value = 2*_sp.sqrt(areas/_sp.pi)
    return value
