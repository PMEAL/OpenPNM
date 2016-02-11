r"""
===============================================================================
Submodule -- throat_diameter
===============================================================================

"""
from . import misc
import scipy as _sp


def weibull(geometry, shape, scale, loc, throat_seed='throat.seed', **kwargs):
    return misc.weibull(geometry=geometry, shape=shape, scale=scale, loc=loc,
                        seeds=throat_seed)
weibull.__doc__ = misc.weibull.__doc__


def normal(geometry, scale, loc, throat_seed='throat.seed', **kwargs):
    return misc.normal(geometry=geometry, scale=scale, loc=loc,
                       seeds=throat_seed)
normal.__doc__ = misc.normal.__doc__


def generic(geometry, func, throat_seed='throat.seed', **kwargs):
    return misc.generic(geometry=geometry, func=func, seeds=throat_seed)
generic.__doc__ = misc.generic.__doc__


def cylinder(geometry, tsd_name, tsd_shape, tsd_loc, tsd_scale,
             throat_seed='throat.seed', tsd_offset=0, **kwargs):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    """
    import scipy.stats as spst
    prob_fn = getattr(spst, tsd_name)
    P = prob_fn(tsd_shape, loc=tsd_loc, scale=tsd_scale)
    value = P.ppf(geometry[throat_seed]) + tsd_offset
    return value


def equivalent_circle(geometry, throat_area='throat.area', **kwargs):
    r"""
    Equivalent circular diameter from throat area
    """
    areas = geometry[throat_area]
    value = 2*_sp.sqrt(areas/_sp.pi)  # 64 bit sqrt doesn't work!
    return value


def minpore(network, geometry, factor=1, **kwargs):
    r"""
    Assign the throat diameter to be equal to the smallest connecting pore
    diameter. If zero (in case of boundaries) take it to be the maximum of
    the connecting pore diameters

    Parameters
    ----------
    factor : float < 1
        A factor between 0 and 1 to further constrict the throat size
        calculcated by the function.

    """
    gTs = geometry.throats()
    nTs = geometry.map_throats(network, gTs)
    pDs = network["pore.diameter"][network["throat.conns"][nTs]]
    value = _sp.amin(pDs, axis=1)*factor
    value[value == 0.0] = _sp.amax(pDs, axis=1)[value == 0.0]
    return value
