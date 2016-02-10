r"""
===============================================================================
Submodule -- Miscillaneous functions
===============================================================================

"""
import scipy as _sp


def random(geometry, element, seed=None, num_range=[0, 1], **kwargs):
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
    value = _sp.random.rand(geometry._count(element),)
    value = value*range_size + range_min
    return value


def weibull(geometry, shape, scale, loc, seeds, **kwargs):
    r"""

    Examples
    --------
    >>> plt.hist(spst.weibull_min.ppf(q=sp.rand(10000), c=2, scale=.0001,
    ...                               loc=0), bins=50)
    """
    value = _sp.stats.weibull_min.ppf(q=geometry[seeds], c=shape, scale=scale,
                                      loc=loc)
    return value


def normal(geometry, scale, loc, seeds, **kwargs):
    r"""

    Examples
    --------
    >>> plt.hist(spst.norm.ppf(q=sp.rand(10000), scale=.0001, loc=0.001),
    ...          bins=50)

    """
    value = _sp.stats.norm.ppf(q=geometry[seeds], scale=scale, loc=loc)
    return value


def generic(geometry, func, seeds, **kwargs):
    r"""
    Accepts an 'rv_frozen' object from the Scipy.stats submodule and returns
    values from the distribution for the given seeds using the ``ppf`` method.

    Examples
    --------
    The folling code illustrates the process of obtaining a 'frozen' Scipy
    stats object, and visualizes the corresponding distribution using a
    Matplotlib histogram:

    .. code-block:: python

        import scipy
        func = scipy.stats.weibull_min(c=2, scale=.0001, loc=0)
        import matplotlib.pylot as plt
        plt.hist(func.ppf(q=scipy.rand(1000), bins=50))

    """
    value = func.ppf(geometry[seeds])
    return value
