r"""

.. autofunction:: openpnm.models.misc.statistical_distributions.random
.. autofunction:: openpnm.models.misc.statistical_distributions.weibull
.. autofunction:: openpnm.models.misc.statistical_distributions.normal
.. autofunction:: openpnm.models.misc.statistical_distributions.generic_distribution

"""
import numpy as np
from scipy import stats
from openpnm.utils import logging
logger = logging.getLogger(__name__)


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

    Returns
    -------
    values : NumPy ndarray
        Array containing uniformly-distributed random numbers.

    """
    range_size = num_range[1] - num_range[0]
    range_min = num_range[0]
    if seed is not None:
        np.random.seed(seed)
    value = np.random.rand(target._count(element),)
    value = value*range_size + range_min
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

    Returns
    -------
    values : NumPy ndarray
        Array containing random numbers based on Weibull distribution.

    Examples
    --------
    The following code illustrates the inner workings of this function,
    which uses the 'weibull_min' method of the scipy.stats module.  This can
    be used to find suitable values of 'shape', 'scale'` and 'loc'.  Note that
    'shape' is represented by 'c' in the actual function call.

    >>> import scipy
    >>> import numpy
    >>> func = scipy.stats.weibull_min(c=1.5, scale=0.0001, loc=0)
    >>> import matplotlib.pyplot as plt
    >>> fig = plt.hist(func.ppf(q=numpy.random.rand(10000)), bins=50)

    """
    import scipy.stats as spts

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

    Returns
    -------
    values : NumPy ndarray
        Array containing normally distributed random numbers.

    Examples
    --------
    The following code illustrates the inner workings of this function,
    which uses the 'norm' method of the scipy.stats module.  This can
    be used to find suitable values of 'scale' and 'loc'.

    >>> import scipy
    >>> import numpy
    >>> func = scipy.stats.norm(scale=.0001, loc=0.001)
    >>> import matplotlib.pyplot as plt
    >>> fig = plt.hist(func.ppf(q=numpy.random.rand(10000)), bins=50)

    """
    import scipy.stats as spts

    seeds = target[seeds]
    value = spts.norm.ppf(q=seeds, scale=scale, loc=loc)
    return value


def generic_distribution(target, seeds, func, **kwargs):
    r"""
    Given the name of a ``scipy.stats`` distribution, uses the ``ppf`` method
    to obtain values.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    seeds : string, optional
        The dictionary key containing random seed values (between 0 and 1)
        to use in the statistical distribution.

    func : str
        The name of the desired ``scipy.stats`` function to use

    keyword arguments
        All additional arguments a passed directly to the ``func``.

    Returns
    -------
    values : NumPy ndarray
        Array containing random numbers based on given ``ppf``

    Notes
    -----
    The ``ppf`` is the reverse lookup of the cumulative density function so
    specifying the y-axis (i.e. seeds between 0 and 1) allows for the reverse
    lookup of the sizes from the x-axis.

    Examples
    --------
    The following code illustrates the process of obtaining a 'frozen' Scipy
    stats object and adding it as a model:

    >>> import scipy
    >>> import numpy
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[3, 3, 3])
    >>> geo = op.geometry.GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)
    >>> geo.add_model(propname='pore.seed',
    ...               model=op.models.geometry.pore_seed.random)

    Now retrieve the stats distribution and add to ``geo`` as a model:

    >>> geo.add_model(propname='pore.size',
    ...               model=op.models.geometry.pore_size.generic_distribution,
    ...               seeds='pore.seed',
    ...               func='weibull_min',
    ...               c=2, scale=0.0001, loc=0)


    >>> import matplotlib.pyplot as plt
    >>> fig = plt.hist(stats_obj.ppf(q=numpy.random.rand(1000)), bins=50)

    """
    func = getattr(stats, func)
    seeds = target[seeds]
    rv_frozen = func(**kwargs)
    return rv_frozen.ppf(seeds)
