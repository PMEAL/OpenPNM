r"""
Pore seed models are use to calculate random numbers for each pore, which can
subsequently be used in statistical distributions to calculate actual pore
sizes.
"""
import numpy as _np
import scipy as _sp
from openpnm.models import misc as _misc
from openpnm.utils import logging as _logging
_logger = _logging.getLogger(__name__)


def random(target, seed=None, num_range=[0, 1]):
    return _misc.random(target, element='pore', seed=seed, num_range=num_range)


random.__doc__ = _misc.random.__doc__


def spatially_correlated(target, weights=None, strel=None):
    r"""
    Generates pore seeds that are spatailly correlated with their neighbors.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

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

    ::

        strel = np.array([[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
                          [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
                          [[0, 0, 0], [0, 0, 0], [0, 0, 0]]])

    Returns
    -------
    values : NumPy ndarray
        Array containing pore seed values.

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
    .. [2] J. Gostick et al, Pore network modeling of fibrous gas diffusion
           layers for polymer electrolyte membrane fuel cells. J Power Sources
           v173, pp277–290 (2007)

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[10, 10, 10])
    >>> Ps, Ts = pn.Ps, pn.Ts
    >>> geom = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)
    >>> mod = op.models.geometry.pore_seed.spatially_correlated
    >>> geom.add_model(propname='pore.seed', model=mod, weights=[2, 2, 2])

    """
    import scipy.ndimage as spim
    network = target.project.network
    # The following will only work on Cubic networks
    x, y, z = network.settings['shape']
    im = _np.random.rand(x, y, z)
    if strel is None:  # Then generate a strel
        if sum(weights) == 0:
            # If weights of 0 are sent, then skip everything and return rands.
            return im.flatten()
        w = _np.array(weights)
        strel = _np.zeros(w*2+1)
        strel[:, w[1], w[2]] = 1
        strel[w[0], :, w[2]] = 1
        strel[w[0], w[1], :] = 1
    im = spim.convolve(im, strel)
    # Convolution is no longer randomly distributed, so fit a gaussian
    # and find it's seeds
    im = (im - _np.mean(im))/_np.std(im)
    im = 1/2*_sp.special.erfc(-im/_np.sqrt(2))
    values = im.flatten()
    values = values[network.pores(target.name)]
    return values
