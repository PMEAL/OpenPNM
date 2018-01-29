r"""
===============================================================================
pore_diameter
===============================================================================

"""
from openpnm.core import logging
from . import misc as _misc
import scipy as _sp
_logger = logging.getLogger(__name__)


def weibull(target, shape, scale, loc, seeds='pore.seed'):
    return _misc.weibull(target=target, shape=shape, scale=scale, loc=loc,
                         seeds=seeds)


weibull.__doc__ = _misc.weibull.__doc__


def normal(target, scale, loc, seeds='pore.seed'):
    return _misc.normal(target=target, scale=scale, loc=loc,
                        seeds=seeds)


normal.__doc__ = _misc.normal.__doc__


def generic(target, func, seeds='pore.seed'):
    return _misc.generic(target=target, func=func, seeds=seeds)


generic.__doc__ = _misc.generic.__doc__


def random(target, seed=None, num_range=[0, 1]):
    r"""
    Assign pore sizes from a random distribution

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
        pore sizes directly; for instance pores between 10 and 100 um can be
        generated with ``num_range = [0.00001, 0.0001]``.
    """
    return _misc.random(target, element='pore', seed=seed, num_range=num_range)


def largest_sphere(target, pore_diameter='pore.diameter', iters=10):
    r"""
    Finds the maximum diameter pore that can be place in each location that
    does not overlap with any neighbors.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    iters : integer
        The number of iterations to perform when searching for maximum
        diameter.  This function iteratively grows pores until they touch
        their nearest neighbor, which is also growing, so this parameter limits
        the maximum number of iterations.  The default is 10, but 5 is usally
        enough.

    Notes
    -----
    This model looks into all pores in the network when finding the diameter.
    This means that when multiple Geometry objects are defined, it will
    consider the diameter of pores on adjacent Geometries. If no diameters
    have been assigned to these neighboring pores it will assume 0.  If
    diameter value are assigned to the neighboring pores AFTER this model is
    run, the pores will overlap.  This can be remedied by running this model
    again.

    """
    network = target.simulation.network
    try:
        D = network[pore_diameter]
        nans = _sp.isnan(D)
        D[nans] = 0.0
    except:
        D = _sp.zeros([network.Np, ], dtype=float)
    Ps = network.pores(target.name)
    C1 = network['pore.coords'][network['throat.conns'][:, 0]]
    C2 = network['pore.coords'][network['throat.conns'][:, 1]]
    L = _sp.sqrt(_sp.sum((C1 - C2)**2, axis=1))
    while iters >= 0:
        iters -= 1
        Lt = L - _sp.sum(D[network['throat.conns']], axis=1)/2
        am = network.create_adjacency_matrix(weights=Lt, fmt='lil')
        D[Ps] = D[Ps] + _sp.array([_sp.amin(row) for row in am.data])[Ps]*0.95
    if _sp.any(D < 0):
        _logger.warning('Negative pore diameters found!  Neighboring pores' +
                        ' must be larger than the pore spacing.')
    return D[network.pores(target.name)]


def equivalent_sphere(target, pore_volume='pore.volume'):
    r"""
    Calculate pore diameter as the diameter of a sphere with an equivalent
    volume.

    Parameters
    ----------
    target : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_volume : string
        The dictionary key containing the pore volume values
    """
    from scipy.special import cbrt
    pore_vols = target[pore_volume]
    value = cbrt(6*pore_vols/_sp.pi)
    return value
