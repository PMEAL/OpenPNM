r"""
===============================================================================
Submodule -- throat_diameter
===============================================================================

"""
import scipy as _sp
import numpy as np


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


def minpore(network, geometry, factor=None, **kwargs):
    r"""
    Assign the throat diameter to be equal to the smallest connecting pore
    diameter. If zero (in case of boundaries) take it to be the maximum of
    the connecting pore diameters
    """
    gTs = geometry.throats()
    nTs = geometry.map_throats(network, gTs)
    pDs = network["pore.diameter"][network["throat.conns"][nTs]]
    value = np.min(pDs, axis=1)
    value[value == 0.0] = np.max(pDs, axis=1)[value == 0.0]
    if factor is not None:
        value *= factor
    return value

def meanpore(network, geometry, factor=None, **kwargs):
    r"""
    Assign the throat diameter to be equal to the mean of the connecting pore
    diameters. If zero (in case of boundaries) take it to be the maximum of
    the connecting pore diameters
    """
    gTs = geometry.throats()
    nTs = geometry.map_throats(network, gTs)
    pDs = network["pore.diameter"][network["throat.conns"][nTs]]
    value = np.mean(pDs, axis=1)
    value[value == 0.0] = np.max(pDs, axis=1)[value == 0.0]
    if factor is not None:
        value *= factor
    return value
