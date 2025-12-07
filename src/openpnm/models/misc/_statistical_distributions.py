import logging
import numpy as np
import scipy.stats as spts
from openpnm.models import _doctxt

logger = logging.getLogger(__name__)


__all__ = [
    'random',
    'weibull',
    'normal',
    'generic_distribution',
    'match_histogram',
]


@_doctxt
def weibull(network, seeds, shape, scale, loc):
    r"""
    Produces values from a Weibull distribution given a set of random numbers.

    Parameters
    ----------
    %(network)s
    seeds : str (dict key)
        %(dict_blurb) seed
    shape : float
        Controls the width or skewness of the distribution. For more
        information on the effect of this parameter refer to the
        corresponding `scipy.stats function
        <https://docs.scipy.org/doc/scipy/reference/stats.html>`_.
    scale : float
        Controls the width of the distribution. For more information on
        the effect of this parameter refer to the corresponding
        scipy.stats function.
    loc : float
            Specifies the central value of ???

    Returns
    -------
    values : ndndarray
        A numpy ndarray containing values following the distribution

    Examples
    --------
    The following code illustrates the inner workings of this function,
    which uses the 'weibull_min' method of the scipy.stats module.  This can
    be used to find suitable values of 'shape', 'scale'` and 'loc'.  Note that
    'shape' is represented by 'c' in the actual function call.

    .. plot::

       import numpy
       import scipy.stats
       import matplotlib.pyplot as plt

       func = scipy.stats.weibull_min(c=1.5, scale=0.0001, loc=0)
       plt.hist(func.ppf(q=numpy.random.rand(10000)), bins=50)

       plt.show()

    """
    seeds = network[seeds]
    value = spts.weibull_min.ppf(q=seeds, c=shape, scale=scale, loc=loc)
    return value


@_doctxt
def normal(network, seeds, mean=None, stddev=None, scale=None, loc=None):
    r"""
    Produces values from a Weibull distribution given a set of random numbers.

    Parameters
    ----------
    %(network)s
    seeds : str (dict key)
        %(dict_blurb) seed
    mean : float
        The mean value of the distribution.  This is referred to as the
        ``loc`` in the scipy.stats function, and this key word is also
        accepted.
    stddev : float
        The standard deviation of the distribution.  This is referred to as
        the  ``scale`` in the scipy.stats function, and this key word is also
        accepted.

    Returns
    -------

    Examples
    --------
    The following code illustrates the inner workings of this function,
    which uses the 'norm' method of the scipy.stats module.  This can
    be used to find suitable values of 'scale' and 'loc'.

    .. plot::

       import numpy
       import scipy.stats
       import matplotlib.pyplot as plt

       func = scipy.stats.norm(scale=.0001, loc=0.001)
       fig = plt.hist(func.ppf(q=numpy.random.rand(10000)), bins=50)

       plt.show()

    """
    scale = stddev if stddev is not None else scale
    loc = mean if mean is not None else loc
    seeds = network[seeds]
    value = spts.norm.ppf(q=seeds, scale=scale, loc=loc)
    return value


@_doctxt
def generic_distribution(network, seeds, func, **kwargs):
    r"""
    Accepts an object from the Scipy.stats submodule and returns
    values from the distribution for the given seeds



    Parameters
    ----------
    %(network)s
    seeds : str (dict key)
        %(dict_blurb) seed
    func : object
        An object from the scipy.stats library. Can be a 'frozen' object, where all
        the parameters were specified upon creation, or a handle to an unintialized
        object. In the latter case the parameters for the distribution should be
        provided as keyword arguments.

    Returns
    -------
    values : ndarray
        An ndarray of either Np or Nt length, with values taken from the supplied
        distribution.

    Examples
    --------
    The following code illustrates the process of obtaining a 'frozen' Scipy
    stats object and adding it as a model:

    .. plot::

       import numpy
       import scipy.stats
       import openpnm as op
       import matplotlib.pyplot as plt

       pn = op.network.Cubic(shape=[3, 3, 3])
       pn.add_model(propname='pore.seed',
                    model=op.models.geometry.pore_seed.random)

       # Now retrieve the stats distribution and add to ``geo`` as a model
       stats_obj = scipy.stats.weibull_min(c=2, scale=.0001, loc=0)
       pn.add_model(propname='pore.size',
                    model=op.models.geometry.pore_size.generic_distribution,
                    seeds='pore.seed',
                    func=stats_obj)

       plt.hist(stats_obj.ppf(q=numpy.random.rand(1000)), bins=50)

       plt.show()

    """
    if hasattr(func, 'freeze'):
        func = func.freeze(**kwargs)
    seeds = network[seeds]
    value = func.ppf(seeds)
    return value


@_doctxt
def random(network, element, seed=None, num_range=[0, 1]):
    r"""
    Create an array of random numbers of a specified size.

    Parameters
    ----------
    %(network)s
    seed : int
        The starting seed value to sent to numpy's random number generator.
        A value of ``None`` means a different distribution is returned each
        time the model is (re)run.
    num_range : list
        A two element list indicating the low and high end of the returned
        numbers.  The default is ``[0, 1]``, but a value of ``[0.1, 0.9]``
        may be useful if these values are to be used in subsequent
        distributions to prevent generating extreme values in the tails.

    Returns
    -------

    """
    range_size = num_range[1] - num_range[0]
    range_min = num_range[0]
    if seed is not None:
        np.random.seed(seed)
    value = np.random.rand(network._count(element),)
    value = value*range_size + range_min
    return value


@_doctxt
def match_histogram(network, bin_centers, bin_heights, element='pore'):
    r"""
    Generate values corresponding to a given histogram

    Parameters
    ----------
    %(network)s
    bin_centers : array_like
        The x-axis of the histogram, such as pore sizes.
    bin_heights : array_like
        The y-axis of the histogram, such as the number of pores of each size.
    element : str
        Controls how many values to generate. Can either be 'pore' or 'throat'.

    Returns
    -------
    values : ndarray
        Values corresponding to ``bin_centers`` generated in proportion to the
        respective ``bin_heights``.

    """
    N = network._count(element)
    h = np.cumsum(bin_heights)
    b = np.digitize(np.random.rand(N)*np.amax(h), bins=h)
    vals = np.array(bin_centers)[b]
    return vals
