r"""
===============================================================================
pore_misc -- miscillaneous and generic functions to apply to pores
===============================================================================

"""
from . import misc as _misc
import scipy as _sp


def constant(target, value):
    r"""
    Assign specified constant value to all pores on the specifed Geometry
    object.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    value : scalar
        The scalar value to be applied

    Notes
    -----
    This function is redundant and could be accomplished with
    geometry['pore.prop'] = value.
    """
    if _sp.shape(value):
        raise Exception('Only scalar values can be assigned')
    value = _sp.ones(target.num_pores(),)*value
    return value


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
    values = _misc.random(target, element='pore', seed=seed,
                          num_range=num_range)
    return values


def neighbor(target, throat_prop='throat.seed', mode='min'):
    r"""
    Adopt a value from the values found in neighboring throats

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_prop : string
        The dictionary key of the array containing the throat property to be
        used in the calculation.  The default is 'throat.seed'.

    mode : string
        Controls how the pore property is calculated.  Options are 'min',
        'max' and 'mean'.
    """
    network = target.project.network
    Ps = target.pores()
    data = target[throat_prop]
    neighborTs = network.find_neighbor_throats(pores=Ps,
                                               flatten=False,
                                               mode='intersection')
    values = _sp.ones((_sp.shape(Ps)[0],))*_sp.nan
    if mode == 'min':
        for pore in Ps:
            values[pore] = _sp.amin(data[neighborTs[pore]])
    if mode == 'max':
        for pore in Ps:
            values[pore] = _sp.amax(data[neighborTs[pore]])
    if mode == 'mean':
        for pore in Ps:
            values[pore] = _sp.mean(data[neighborTs[pore]])
    return values
