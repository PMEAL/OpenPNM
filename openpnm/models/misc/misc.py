import numpy as np
import scipy.stats as spts


def constant(target, value):
    r"""
    Places a constant value into the target object

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    value : scalar
        The numerical value to apply

    Notes
    -----
    This model is mostly useless and for testing purposes, but might be used
    to 'reset' an array back to a default value.
    """
    return value


def product(target, props):
    r"""
    Calculates the product of multiple property values

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    props : list of strings
        The names of each property to be multiplied

    """
    value = 1.0
    for item in props:
        value *= target[item]
    return value


def random(target, element, seed=None, num_range=[0, 1]):
    r"""
    Create an array of random numbers of a specified size.

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
        numbers.

    """
    range_size = num_range[1] - num_range[0]
    range_min = num_range[0]
    np.random.seed(seed)
    value = np.random.rand(target._count(element),)
    value = value*range_size + range_min
    return value


def scaled(target, prop, factor):
    r"""
    Scales an existing value by a factor.  Useful for constricting some throat
    property.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    prop : string
        The dictionary key of the array containing the values to be scaled.

    factor : scalar
        The factor by which the values should be scaled.
    """
    value = target[prop]*factor
    return value


def linear(target, m, b, prop):
    r"""
    Calculates a property as a linear function of a given property

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    m, b : floats
        Slope and intercept of the linear corelation

    prop : string
        The dictionary key containing the independent variable or phase
        property to be used in the correlation.

    """
    x = target[prop]
    value = m*x + b
    return value


def polynomial(target, a, prop, **kwargs):
    r"""
    Calculates a property as a polynomial function of a given property

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    a : array_like
        A list containing the polynomial coefficients, where element 0 in the
        list corresponds to a0 and so on.  Note that no entries can be skipped
        so 0 coefficients must be sent as 0.

    poreprop : string
        The dictionary key containing the independent variable or phase
        property to be used in the correlation.

    """
    x = target[prop]
    value = 0.0
    for i in range(0, len(a)):
        value += a[i]*x**i
    return value


def weibull(target, seeds, shape, scale, loc):
    r"""
    Produces values from a Weibull distribution given a set of random numbers.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    seeds : string, optional
        The dictionary key on the Geometry object containing random seed values
        (between 0 and 1) to use in the statistical distribution.

    shape : float
        This controls the skewness of the distribution, with 'shape' < 1 giving
        values clustered on the low end of the range with a long tail, and
        'shape' > 1 giving a more symmetrical distribution.

    scale : float
        This controls the width of the distribution with most of values falling
        below this number.

    loc : float
        Applies an offset to the distribution such that the smallest values are
        above this number.

    Examples
    --------
    The following code illustrates the inner workings of this function,
    which uses the 'weibull_min' method of the scipy.stats module.  This can
    be used to find suitable values of 'shape', 'scale'` and 'loc'.  Note that
    'shape' is represented by 'c' in the actual function call.

    >>> import scipy
    >>> func = scipy.stats.weibull_min(c=1.5, scale=0.0001, loc=0)
    >>> import matplotlib.pyplot as plt
    >>> fig = plt.hist(func.ppf(q=scipy.rand(10000)), bins=50)

    """
    seeds = target[seeds]
    value = spts.weibull_min.ppf(q=seeds, c=shape, scale=scale, loc=loc)
    return value


def normal(target, seeds, scale, loc):
    r"""
    Produces values from a Weibull distribution given a set of random numbers.

    Parameters
    ----------
    target : OpenPNM Object
        The object with which this function as associated.  This argument
        is required to (1) set number of values to generate (geom.Np or
        geom.Nt) and (2) provide access to other necessary values
        (i.e. geom['pore.seed']).

    seeds : string, optional
        The dictionary key on the Geometry object containing random seed values
        (between 0 and 1) to use in the statistical distribution.

    scale : float
        The standard deviation of the Normal distribution

    loc : float
        The mean of the Normal distribution

    Examples
    --------
    The following code illustrates the inner workings of this function,
    which uses the 'norm' method of the scipy.stats module.  This can
    be used to find suitable values of 'scale' and 'loc'.

    >>> import scipy
    >>> func = scipy.stats.norm(scale=.0001, loc=0.001)
    >>> import matplotlib.pyplot as plt
    >>> fig = plt.hist(func.ppf(q=scipy.rand(10000)), bins=50)

    """
    seeds = target[seeds]
    value = spts.norm.ppf(q=seeds, scale=scale, loc=loc)
    return value


def generic(target, seeds, func):
    r"""
    Accepts an 'rv_frozen' object from the Scipy.stats submodule and returns
    values from the distribution for the given seeds using the ``ppf`` method.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    seeds : string, optional
        The dictionary key on the Geometry object containing random seed values
        (between 0 and 1) to use in the statistical distribution.

    func : object
        An 'rv_frozen' object from the Scipy.stats library with all of the
        parameters pre-specified.

    Examples
    --------
    The following code illustrates the process of obtaining a 'frozen' Scipy
    stats object, and visualizes the corresponding distribution using a
    Matplotlib histogram:

    >>> import scipy
    >>> func = scipy.stats.weibull_min(c=2, scale=.0001, loc=0)
    >>> import matplotlib.pyplot as plt
    >>> fig = plt.hist(func.ppf(q=scipy.rand(1000)), bins=50)

    """
    seeds = target[seeds]
    value = func.ppf(seeds)
    return value


def from_neighbor_throats(target, throat_prop='throat.seed', mode='min'):
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
    values = np.ones((np.shape(Ps)[0],))*np.nan
    if mode == 'min':
        for pore in Ps:
            values[pore] = np.amin(data[neighborTs[pore]])
    if mode == 'max':
        for pore in Ps:
            values[pore] = np.amax(data[neighborTs[pore]])
    if mode == 'mean':
        for pore in Ps:
            values[pore] = np.mean(data[neighborTs[pore]])
    return values


def from_neighbor_pores(target, pore_prop='pore.seed', mode='min'):
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
        value = np.amin(pvalues, axis=1)
    if mode == 'max':
        value = np.amax(pvalues, axis=1)
    if mode == 'mean':
        value = np.mean(pvalues, axis=1)
    return value
