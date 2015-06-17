r"""
===============================================================================
Submodule -- throat_seeds
===============================================================================

"""
import scipy as _sp


def random(geometry, seed=None, num_range=[0, 1], **kwargs):
    r"""
    Assign random number to throats, for use in statistical distributions that
    return pore size

    Parameters
    ----------
    seed : int
        The starting seed value to send to Scipy's random number generator.
        The default is None, which means different distribution is returned
        each time the model is run.

    num_range : list
        A two element list indicating the low and high end of the returned
        numbers.
    """
    range_size = num_range[1] - num_range[0]
    range_min = num_range[0]
    _sp.random.seed(seed)
    value = _sp.random.rand(geometry.num_throats(),)
    value = value*range_size + range_min
    return value


def neighbor(geometry, network, pore_prop='pore.seed', mode='min', **kwargs):
    r"""
    Adopt a value based on the values in the neighboring pores

    Parameters
    ----------
    mode : string
        Indicates how to select the values from the neighboring pores.  The
        options are:

        - min : (Default) Uses the minimum of the value found in the neighbors
        - max : Uses the maximum of the values found in the neighbors
        - mean : Uses an average of the neighbor values

    pore_prop : string
        The dictionary key containing the pore property to be used.
    """
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
