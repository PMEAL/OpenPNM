r"""
===============================================================================
throat_misc --  Miscillaneous and generic functions to apply to throats
===============================================================================

"""
import scipy as _sp
from . import misc as _misc


def random(target, seed=None, num_range=[0, 1]):
    r"""
    Generate an array of random numbers using the Scipy random number
    generator.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    seed : int
        The seed value to initialize the random number generator.  The default
        is None which means different numbers will be produced each time this
        function is run.

    num_range : list
        The min and max values of the returned random numbers.  The values will
        be uniformly distributed over this range.

    Notes
    -----
    This method seems trivial but it can be useful to add this as a pore-scale
    model so that random number can be regenerated upon demand to create
    new realization of the pore and throat size distributions (i.e. assuming
    ``seed`` is not set to None).
    """
    values = _misc.random(target, element='throat', seed=seed,
                          num_range=num_range)
    return values


def neighbor(target, pore_prop='pore.seed', mode='min'):
    r"""
    Adopt a value based on the values in neighboring pores

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_prop : string
        The dictionary key to the array containing the pore property to be
        used in the calculation.  Default is 'pore.seed'.

    mode : string
        Controls how the throat property is calculated.  Options are 'min',
        'max' and 'mean'.

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
