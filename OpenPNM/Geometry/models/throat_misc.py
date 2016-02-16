r"""
===============================================================================
throat_misc --  Miscillaneous and generic functions to apply to throats
===============================================================================

"""
import scipy as _sp


def random(geometry, seed=None, num_range=[0, 1], **kwargs):
    r"""
    Assign random number to throats


    """
    range_size = num_range[1] - num_range[0]
    range_min = num_range[0]
    _sp.random.seed(seed=seed)
    value = _sp.random.rand(geometry.num_throats(),)
    value = value*range_size + range_min
    return value


def neighbor(geometry, pore_prop='pore.seed', mode='min', **kwargs):
    r"""
    Adopt a value based on the values in neighboring pores

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The object containing the ``pore_prop`` to be used.

    pore_prop : string
        The dictionary key to the array containing the pore property to be
        used in the calculation.  Default is 'pore.seed'.

    mode : string
        Controls how the throat property is calculated.  Options are 'min',
        'max' and 'mean'.

    """
    network = geometry._net
    throats = network.throats(geometry.name)
    P12 = network.find_connected_pores(throats)
    pvalues = network[pore_prop][P12]
    if mode == 'min':
        value = _sp.amin(pvalues, axis=1)
    if mode == 'max':
        value = _sp.amax(pvalues, axis=1)
    if mode == 'mean':
        value = _sp.mean(pvalues, axis=1)
    return value
