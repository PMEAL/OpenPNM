# -*- coding: utf-8 -*-
r"""
===============================================================================
pore_seed -- Methods for generating fields of values for use as seeds in
statistical pore size distributions
===============================================================================

"""
import scipy as _sp


def random(geometry, seed=None, num_range=[0, 1], **kwargs):
    r"""
    Assign random number to pores, for use in statistical distributions that
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
    range_size = num_range[1]-num_range[0]
    range_min = num_range[0]
    _sp.random.seed(seed=seed)
    value = _sp.random.rand(geometry.num_pores(),)
    value = value*range_size + range_min
    return value


def distance_from_inclusion(geometry, p, **kwargs):
    r"""
    Genrate spatially correlated pore seeds by calculating distance from random
    locations (inclusions) in the domain

    Parameters
    ----------
    p : float
        The fraction of pores in the domain that are set as 'seeds' for the
        distance calculation

    Returns
    -------
    A list of distance values (in voxels) between each pore and it nearest
    seed pore.  A list of voxel distances is returned rather than normalized
    seeds between 0:1 so that the user can manipulate the map as desired, by
    applying desired thresholds and/or scaling to get 0:1 seeds.

    Notes
    -----
    - This method uses image analysis tools, so only works on Cubic networks
    - At present the result contains edge artifacts since no inclusions are present
      beyond the image boundary

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.Cubic(shape=[50,50,50])
    >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,pores=pn.Ps,throats=pn.Ts)
    >>> model = OpenPNM.Geometry.models.pore_seed.distance_from_inclusion
    >>> geom.add_model(propname='pore.seed', model=model, p=0.001)
    >>> im = pn.asarray(geom['pore.seed'])

    Visualizing the end result can be done with:

    .. code-block:: python

        matplotlib.pyplot.imshow(im[:,25,:],interpolation='none')

    """
    import scipy.ndimage as _spim
    net = geometry._net
    # The following will only work on Cubic networks
    x = net._shape[0]
    y = net._shape[1]
    z = net._shape[2]
    img = _sp.rand(x, y, z) > p
    # Pad image by tiling
    a = _sp.tile(img, [3, 3, 3])
    b = a[x:-x, y:-y, z:-z]
    # Perform distance transform
    img = _spim.distance_transform_bf(b)
    # Convert back to pore-list
    values = img.flatten()
    values = values[geometry.map_pores(target=net, pores=geometry.Ps)]
    return values


def spatially_correlated(geometry, network, weights=None, strel=None, **kwargs):
    r"""
    Generates pore seeds that are spatailly correlated with their neighbors.

    Parameters
    ----------
    weights : list of ints, optional
        The [Nx,Ny,Nz] distances (in number of pores) in each direction that
        should be correlated.

    strel : array_like, optional (in place of weights)
        The option allows full control over the spatial correlation pattern by
        specifying the structuring element to be used in the convolution.

        The array should be a 3D array containing the strength of correlations
        in each direction.  Nonzero values indicate the strength, direction
        and extent of correlations.  The following would achieve a basic
        correlation in the z-direction:

        strel = sp.array([[[0, 0, 0], [0, 0, 0], [0, 0, 0]],\
                            [[0, 0, 0], [1, 1, 1], [0, 0, 0]],\
                            [[0, 0, 0], [0, 0, 0], [0, 0, 0]]])

    Notes
    -----
    This approach uses image convolution to replace each pore seed in the
    geoemtry with a weighted average of those around it.  It then converts the
    new seeds back to a random distribution by assuming they new seeds are
    normally distributed.

    Because is uses image analysis tools, it only works on Cubic networks.

    This is the appproached used by Gostick et al [2]_ to create an anistropic
    gas diffusion layer for fuel cell electrodes.

    References
    ----------
    .. [2] J. Gostick et al, Pore network modeling of fibrous gas diffusion layers
           for polymer electrolyte membrane fuel cells. J Power Sources 173 (2007)
           277–290

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.Cubic(shape=[50, 50, 50])
    >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,
    ...                                         pores=pn.Ps,
    ...                                         throats=pn.Ts)
    >>> geom.add_model(propname='pore.seed',
    ...               model=OpenPNM.Geometry.models.pore_seed.spatially_correlated,
    ...               weights=[2, 2, 2])
    >>> im = pn.asarray(geom['pore.seed'])

    Visualizing the end result can be done with:

    .. code-block:: python

        matplotlib.pyplot.imshow(im[:, 25, :],interpolation='none')

    """
    import scipy.ndimage as spim
    import scipy.stats as spst
    # The following will only work on Cubic networks
    x = network._shape[0]
    y = network._shape[1]
    z = network._shape[2]
    im = _sp.rand(x, y, z)
    if strel is None:  # Then generate a strel
        if sum(weights) == 0:
            # If weights of 0 are sent, then skip everything and return rands.
            return im.flatten()
        w = _sp.array(weights)
        strel = _sp.zeros(w*2+1)
        strel[:, w[1], w[2]] = 1
        strel[w[0], :, w[2]] = 1
        strel[w[0], w[1], :] = 1
    im = spim.convolve(im, strel)
    # Convolution is no longer randomly distributed, so fit a gaussian
    # and find it's seeds
    temp = im.flatten()
    x_mean = _sp.mean(temp)
    x_sigma = _sp.sqrt(1/(temp.size-1)*_sp.sum((temp - x_mean)**2))
    fn1 = spst.norm(loc=x_mean, scale=x_sigma)
    values = fn1.cdf(temp)
    values = values[geometry.map_pores(target=network, pores=geometry.Ps)]
    return values
