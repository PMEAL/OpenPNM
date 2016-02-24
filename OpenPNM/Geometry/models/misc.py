r"""
===============================================================================
Submodule -- Miscillaneous functions
===============================================================================

"""
import scipy as _sp
import scipy.stats as _spts


def random(N, seed=None, num_range=[0, 1], **kwargs):
    r"""
    Create an array of random numbers of a specified size.

    Parameters
    ----------
    N : int
        The length of the random array to be generated.  This should correspond
        to the number of pores or throats on the Geometry object.  For instance
        N = geom.Nt or N = geom.Np.

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
    value = _sp.random.rand(N,)
    value = value*range_size + range_min
    return value


def scaled(geometry, prop, factor, **kwargs):
    r"""
    Scales an existing value by a factor.  Useful for constricting some throat
    property.

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The object with which this function as associated.  This argument
        is required to (1) set number of values to generate (geom.Np or
        geom.Nt) and (2) provide access to other necessary values
        (i.e. geom['pore.seed']).

    prop : string
        The dictionary key of the array containing the values to be scaled.

    factor : scalar
        The factor by which the values should be scaled.
    """
    value = geometry[prop]*factor
    return value


def weibull(geometry, shape, scale, loc, seeds, **kwargs):
    r"""
    Produces values from a Weibull distribution given a set of random numbers.

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The object with which this function as associated.  This argument
        is required to (1) set number of values to generate (geom.Np or
        geom.Nt) and (2) provide access to other necessary values
        (i.e. geom['pore.seed']).

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

    seeds : string, optional
        The dictionary key on the Geometry object containing random seed values
        (between 0 and 1) to use in the statistical distribution.  If none is
        specified, then an array of random numbers will be automatically
        generated and stored on the Geometry object.

    Examples
    --------
    The following code illustrates the inner workings of this function,
    which uses the 'weibull_min' method of the scipy.stats module.  This can
    be used to find suitable values of 'shape', 'scale'` and 'loc'.  Note that
    'shape' is represented by 'c' in the actual function call.

    .. code-block::

        import scipy.stats as spst
        import matplotlib.pyplot as plt
        x = spst.weibull_min.ppf(q=sp.rand(10000), c=1.5, scale=0.0001, loc=0)
        plt.hist(x, bins=50)

    """
    seeds = geometry[seeds]
    value = _spts.weibull_min.ppf(q=seeds, c=shape, scale=scale, loc=loc)
    return value


def normal(geometry, scale, loc, seeds, **kwargs):
    r"""
    Produces values from a Weibull distribution given a set of random numbers.

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The object with which this function as associated.  This argument
        is required to (1) set number of values to generate (geom.Np or
        geom.Nt) and (2) provide access to other necessary values
        (i.e. geom['pore.seed']).

    scale : float
        The standard deviation of the Normal distribution

    loc : float
        The mean of the Normal distribution

    seeds : string, optional
        The dictionary key on the Geometry object containing random seed values
        (between 0 and 1) to use in the statistical distribution.  If none is
        specified, then an array of random numbers will be automatically
        generated and stored on teh Geometry object.

    Examples
    --------
    The following code illustrates the inner workings of this function,
    which uses the 'norm' method of the scipy.stats module.  This can
    be used to find suitable values of 'scale' and 'loc'.

    .. code-block::

        import scipy.stats as spst
        import matplotlib.pyplot as plt
        x = spst.norm.ppf(q=sp.rand(10000), scale=.0001, loc=0.001)
        plt.hist(x, bins=50)

    """
    seeds = geometry[seeds]
    value = _spts.norm.ppf(q=seeds, scale=scale, loc=loc)
    return value


def generic(geometry, func, seeds, **kwargs):
    r"""
    Accepts an 'rv_frozen' object from the Scipy.stats submodule and returns
    values from the distribution for the given seeds using the ``ppf`` method.

    Parameters
    ----------
    func : object
        An 'rv_frozen' object from the Scipy.stats library with all of the
        parameters pre-specified.

    seeds : string, optional
        The dictionary key on the Geometry object containing random seed values
        (between 0 and 1) to use in the statistical distribution.  If none is
        specified, then an array of random numbers will be automatically
        generated and stored on teh Geometry object.

    Examples
    --------
    The following code illustrates the process of obtaining a 'frozen' Scipy
    stats object, and visualizes the corresponding distribution using a
    Matplotlib histogram:

    .. code-block:: python

        import scipy
        func = scipy.stats.weibull_min(c=2, scale=.0001, loc=0)
        import matplotlib.pylot as plt
        plt.hist(func.ppf(q=scipy.rand(1000), bins=50))

    """
    seeds = geometry[seeds]
    value = func.ppf(seeds)
    return value
