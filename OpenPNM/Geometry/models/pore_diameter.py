r"""
===============================================================================
pore_diameter
===============================================================================

"""
from OpenPNM.Base import logging
from . import misc as _misc
import scipy as _sp
_logger = logging.getLogger()


def weibull(geometry, shape, scale, loc, seeds='pore.seed', **kwargs):
    if seeds not in geometry:
        geometry['pore.seed'] = _sp.rand(geometry.Np,)
    return _misc.weibull(geometry=geometry, shape=shape, scale=scale, loc=loc,
                         seeds=seeds)
weibull.__doc__ = _misc.weibull.__doc__


def normal(geometry, scale, loc, seeds='pore.seed', **kwargs):
    if seeds not in geometry:
        geometry['pore.seed'] = _sp.rand(geometry.Np,)
    return _misc.normal(geometry=geometry, scale=scale, loc=loc,
                        seeds=seeds)
normal.__doc__ = _misc.normal.__doc__


def generic(geometry, func, seeds='pore.seed', **kwargs):
    if seeds not in geometry:
        geometry['pore.seed'] = _sp.rand(geometry.Np,)
    return _misc.generic(geometry=geometry, func=func, seeds=seeds)
generic.__doc__ = _misc.generic.__doc__


def random(geometry, seed=None, num_range=[0, 1], **kwargs):
    r"""
    Assign pore sizes from a random distribution

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
        pore sizes directly; for instance pores between 10 and 100 um can be
        generated with ``num_range = [0.00001, 0.0001]``.
    """
    return _misc.random(N=geometry.Np, seed=seed, num_range=num_range)


def largest_sphere(geometry, network, iters=10, **kwargs):
    r"""
    Finds the maximum diameter pore that can be place in each location that
    does not overlap with any neighbors.

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    network : OpenPNM Network Object
        The Netowrk object is required to lookup the connectivty of each pore
        to find the neighbors and subsequently their separation distance.

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
    try:
        D = network['pore.diameter']
        nans = _sp.isnan(D)
        D[nans] = 0.0
    except:
        D = _sp.zeros([network.Np, ], dtype=float)
    Ps = network.pores(geometry.name)
    C1 = network['pore.coords'][network['throat.conns'][:, 0]]
    C2 = network['pore.coords'][network['throat.conns'][:, 1]]
    L = _sp.sqrt(_sp.sum((C1 - C2)**2, axis=1))
    while iters >= 0:
        iters -= 1
        Lt = L - _sp.sum(D[network['throat.conns']], axis=1)/2
        am = network.create_adjacency_matrix(data=Lt, sprsfmt='lil',
                                             dropzeros=False)
        D[Ps] = D[Ps] + _sp.array([_sp.amin(row) for row in am.data])[Ps]*0.95
    if _sp.any(D < 0):
        _logger.warning('Negative pore diameters found!  Neighboring pores' +
                        ' must be larger than the pore spacing.')
    return D[network.pores(geometry.name)]


def sphere(geometry, psd_name, psd_shape, psd_loc, psd_scale,
           pore_seed='pore.seed', psd_offset=0, **kwargs):
    r"""
    Calculate pore diameter from given seed values.

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    psd_name : string
        The name of the statistical distribution to use. This model uses the
        Scipy.stats module, so any of the distributions available there are
        suitable options.

    psd_shape, loc and scale : scalars
        The parameters to send to the selected statistics model.  Most of the
        Scipy.stats models accept the same keyword arguments.  Note that the
        psd_ prefix is added by OpenPNM to indicate 'pore size distribution'.

    psd_offset : scalar
        Controls the minimum value in the pore size distribution by shifting
        the entire set of values by the given offset.  This is useful for
        avoiding pore sizes too close to zero.

    Notes
    -----
    This pore-scale model is deprecated.  Use ``weibull``, ``normal`` or
    ``generic`` to get produce pore sizes distributions.

    """
    import scipy.stats as spst
    prob_fn = getattr(spst, psd_name)
    P = prob_fn(psd_shape, loc=psd_loc, scale=psd_scale)
    value = P.ppf(geometry[pore_seed]) + psd_offset
    return value


def equivalent_sphere(geometry, pore_volume='pore.volume', **kwargs):
    r"""
    Calculate pore diameter as the diameter of a sphere with an equivalent
    volume.

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_volume : string
        The dictionary key containing the pore volume values
    """
    from scipy.special import cbrt
    pore_vols = geometry[pore_volume]
    value = cbrt(6*pore_vols/_sp.pi)
    return value


def equivalent_cube(geometry, pore_volume='pore.volume', **kwargs):
    r"""
    Calculate pore diameter as the width of a cube with an equivalent volume.

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_volume : string
        The dictionary key containing the pore volume values
    """
    from scipy.special import cbrt
    pore_vols = geometry[pore_volume]
    value = cbrt(pore_vols)
    return value


def centroids(geometry, throat_centroid='throat.centroid',
              pore_centroid='pore.centroid', **kwargs):
    r"""
    Calculate the diameter representing an inclosed sphere. The maximum is very
    difficult to caluclate for irregular polygons with more than 4 faces so an
    average distance from the pore centroid to the throat centroid is an
    approximation.

    Parameters
    ----------
    geometry : OpenPNM Geometry object
        The Geometry object with which this model is associated.  This is
        needed to access the pore and throat centroid values.

    pore_centroid and throat_centroid : string
        Dictionary keys to the arrays containing the pore and throat centroid
        coordinates.
    """
    network = geometry._net
    Np = geometry.num_pores()
    value = _sp.zeros(Np)
    pore_map = geometry.map_pores(target=network,
                                  pores=geometry.pores(),
                                  return_mapping=True)

    for i, net_pore in enumerate(pore_map['target']):
        geom_pore = pore_map['source'][i]
        net_throats = geometry._net.find_neighbor_throats(net_pore)
        geom_throats = geometry._net.map_throats(target=geometry,
                                                 throats=net_throats,
                                                 return_mapping=True)['target']
        tcs = geometry[throat_centroid][geom_throats]
        pc = geometry[pore_centroid][geom_pore]
        value[geom_pore] = _sp.mean(_sp.sqrt(((tcs-pc)*(tcs-pc))[:, 0] +
                                             ((tcs-pc)*(tcs-pc))[:, 1] +
                                             ((tcs-pc)*(tcs-pc))[:, 2]))*2
    return value
