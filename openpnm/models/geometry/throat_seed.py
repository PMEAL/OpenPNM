r"""
===============================================================================
Submodule -- throat_seeds
===============================================================================

"""
from . import misc as _misc
import scipy as _sp


def random(target, seed=None, num_range=[0, 1]):
    return _misc.random(target, element='throat', seed=seed,
                        num_range=num_range)


random.__doc__ = _misc.random.__doc__


def neighbor(target, pore_prop='pore.seed', mode='min'):
    r"""
    Adopt a value based on the values in the neighboring pores

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    mode : string
        Indicates how to select the values from the neighboring pores.  The
        options are:

        - min : (Default) Uses the minimum of the value found in the neighbors
        - max : Uses the maximum of the values found in the neighbors
        - mean : Uses an average of the neighbor values

    pore_prop : string
        The dictionary key containing the pore property to be used.
    """
    network = target.project.network
    throats = network.throats(target.name)
    P12 = network.find_connected_pores(throats)
    pvalues = network[pore_prop][P12]
    if mode == 'min':
        value = _sp.amin(pvalues, axis=1)
    if mode == 'max':
        value = _sp.amax(pvalues, axis=1)
    if mode == 'mean':
        value = _sp.mean(pvalues, axis=1)
    return value
