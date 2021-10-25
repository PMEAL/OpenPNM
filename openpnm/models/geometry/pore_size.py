r"""
The pore-size models in this sub-module are used to apply desired pore-size
distributions to a pore network.  Most of the models accept pore seeds, and
then look-up pore-sizes from cumuative distribtions functions.  There is also
a model for finding the largest possible sphere that can be placed on each
site.
"""
from openpnm.utils import logging as _logging
from openpnm.models import misc as _misc
import numpy as _np
_logger = _logging.getLogger(__name__)


def weibull(target, shape, scale, loc, seeds='pore.seed'):
    return _misc.weibull(target=target, shape=shape, scale=scale, loc=loc,
                         seeds=seeds)


weibull.__doc__ = _misc.weibull.__doc__


def normal(target, scale, loc, seeds='pore.seed'):
    return _misc.normal(target=target, scale=scale, loc=loc,
                        seeds=seeds)


normal.__doc__ = _misc.normal.__doc__


def random(target, seed=None, num_range=[0, 1]):
    return _misc.random(target=target, element='pore', seed=seed,
                        num_range=num_range)


random.__doc__ = _misc.random.__doc__


def generic_distribution(target, func, seeds='pore.seed'):
    return _misc.generic_distribution(target=target, func=func, seeds=seeds)


generic_distribution.__doc__ = _misc.generic_distribution.__doc__


def from_neighbor_throats(target, prop, mode='max'):
    return _misc.from_neighbor_throats(target=target,
                                       prop=prop,
                                       mode=mode)


from_neighbor_throats.__doc__ = _misc.from_neighbor_throats.__doc__


def largest_sphere(target, fixed_diameter='pore.fixed_diameter', iters=5):
    r"""
    Finds the maximum diameter pore that can be placed in each location without
    overlapping any neighbors.

    This method iteratively expands pores by increasing their diameter to
    encompass half of the distance to the nearest neighbor.  If the neighbor
    is not growing because it's already touching a different neighbor, then
    the given pore will never quite touch this neighbor.  Increating the value
    of ``iters`` will get it closer, but it's case of
    [Zeno's paradox](https://en.wikipedia.org/wiki/Zeno%27s_paradoxes) with
    each step cutting the remaining distance in half

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    fixed_diameter : string
        The dictionary key containing the pore diameter values already
        assigned to network, if any.  If not provided a starting value is
        assumed as half-way to the nearest neighbor.

    iters : integer
        The number of iterations to perform when searching for maximum
        diameter.  This function iteratively grows pores until they touch
        their nearest neighbor, which is also growing, so this parameter limits
        the maximum number of iterations.  The default is 10, but 5 is usally
        enough.

    Returns
    -------
    D : NumPy ndarray
        Array containing pore diameter values.

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
    network = target.project.network
    P12 = network['throat.conns']
    C1 = network['pore.coords'][network['throat.conns'][:, 0]]
    C2 = network['pore.coords'][network['throat.conns'][:, 1]]
    L = _np.sqrt(_np.sum((C1 - C2)**2, axis=1))
    try:
        # Fetch any existing pore diameters on the network
        D = network[fixed_diameter]
        # Set any unassigned values (nans) to 0
        D[_np.isnan(D)] = 0
    except KeyError:
        _logger.info('Pore sizes not present, calculating starting values '
                     + 'as half-way to the nearest neighbor')
        D = _np.inf*_np.ones([network.Np, ], dtype=float)
        _np.minimum.at(D, P12[:, 0], L)
        _np.minimum.at(D, P12[:, 1], L)
    while iters >= 0:
        iters -= 1
        Lt = L - _np.sum(D[P12], axis=1)/2
        Dadd = _np.ones_like(D)*_np.inf
        _np.minimum.at(Dadd, P12[:, 0], Lt)
        _np.minimum.at(Dadd, P12[:, 1], Lt)
        D += Dadd
    if _np.any(D < 0):
        _logger.info('Negative pore diameters found!  Neighboring pores are '
                     + 'larger than the pore spacing.')
    return D[network.pores(target.name)]


def equivalent_diameter(target, pore_volume='pore.volume',
                        pore_shape='sphere'):
    r"""
    Calculates the diameter of a sphere or edge-length of a cube with same
    volume as the pore.

    Parameters
    ----------
    target : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_volume : string
        The dictionary key containing the pore volume values

    pore_shape : string
        The shape of the pore body to assume when back-calculating from
        volume.  Options are 'sphere' (default) or 'cube'.

    Returns
    -------
    D : NumPy ndarray
        Array containing pore diameter values.

    """
    from scipy.special import cbrt
    pore_vols = target[pore_volume]
    if pore_shape.startswith('sph'):
        value = cbrt(6*pore_vols/_np.pi)
    elif pore_shape.startswith('cub'):
        value = cbrt(pore_vols)
    return value
