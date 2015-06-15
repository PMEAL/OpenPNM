r"""
===============================================================================
pore_diameter
===============================================================================

"""
import scipy as _sp
import numpy as np

from scipy.spatial import ConvexHull


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


def centroids(network, geometry, **kwargs):
    r"""
    Calculate the diameter representing an inclosed sphere. The maximum is very
    difficult to caluclate for irregular polygons with more than 4 faces so an
    average distance from the pore centroid to the throat centroid is an
    approximation
    """
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
                                                 return_mapping=False)
        tcs = geometry["throat.centroid"][geom_throats]
        pc = geometry["pore.centroid"][geom_pore]
        value[geom_pore] = _sp.mean(_sp.sqrt(((tcs-pc)*(tcs-pc))[:, 0] +
                                             ((tcs-pc)*(tcs-pc))[:, 1] +
                                             ((tcs-pc)*(tcs-pc))[:, 2]))*2
    return value
