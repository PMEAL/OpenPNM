r"""
===============================================================================
Submodule -- throat_diameter
===============================================================================

"""
from . import misc as _misc
import scipy as _sp


def weibull(geometry, shape, scale, loc, seeds=None, **kwargs):
    if seeds not in geometry:
        geometry['throat.seed'] = _sp.rand(geometry.Nt,)
    return _misc.weibull(geometry=geometry, shape=shape, scale=scale, loc=loc,
                         seeds=seeds)
weibull.__doc__ = _misc.weibull.__doc__


def normal(geometry, scale, loc, seeds=None, **kwargs):
    if seeds not in geometry:
        geometry['throat.seed'] = _sp.rand(geometry.Nt,)
    return _misc.normal(geometry=geometry, scale=scale, loc=loc,
                        seeds=seeds)
normal.__doc__ = _misc.normal.__doc__


def generic(geometry, func, seeds=None, **kwargs):
    if seeds not in geometry:
        geometry['throat.seed'] = _sp.rand(geometry.Nt,)
    return _misc.generic(geometry=geometry, func=func, seeds=seeds)
generic.__doc__ = _misc.generic.__doc__


def random(geometry, seed=None, num_range=[0, 1], **kwargs):
    r"""
    Assign throat sizes from a random

    Parameters
    ----------
    geometry : OpenPNM Geometry object
        The Geometry object to which this model will apply.  This is necessary
        to determine the length of the array to generate.

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
    return _misc.random(N=geometry.Nt, seed=seed, num_range=num_range)


def cylinder(geometry, tsd_name, tsd_shape, tsd_loc, tsd_scale,
             throat_seed='throat.seed', tsd_offset=0, **kwargs):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    tsd_name : string
        The name of the statistical distribution to use. This model uses the
        Scipy.stats module, so any of the distributions available there are
        suitable options.

    tsd_shape, loc and scale : scalars
        The parameters to send to the selected statistics model.  Most of the
        Scipy.stats models accept the same keyword arguments.  Note that the
        psd_ prefix is added by OpenPNM to indicate 'pore size distribution'.

    tsd_offset : scalar
        Controls the minimum value in the pore size distribution by shifting
        the entire set of values by the given offset.  This is useful for
        avoiding pore sizes too close to zero.

    Notes
    -----
    This pore-scale model is deprecated.  Use ``weibull``, ``normal`` or
    ``generic`` to get produce pore sizes distributions.
    """
    import scipy.stats as spst
    prob_fn = getattr(spst, tsd_name)
    P = prob_fn(tsd_shape, loc=tsd_loc, scale=tsd_scale)
    value = P.ppf(geometry[throat_seed]) + tsd_offset
    return value


def equivalent_circle(geometry, throat_area='throat.area', **kwargs):
    r"""
    Calculates the diameter of a cirlce with same area as the throat.

    Parameters
    ----------
    geometry : OpenPNM Geometry object
        The Geometry object which contains the throat area values

    thorat_area : string
        The dictionary key to the throat area values
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
