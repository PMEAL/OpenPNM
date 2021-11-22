import numpy as _np


def pore_coords(target):
    r"""
    Calculate throat centroid values by averaging adjacent pore coordinates

    Parameters
    ----------
    target : OpenPNM Base object
        Object with which this model is associated. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.

    Returns
    -------
    value : ndarray
        A numpy ndarray containing throat centroid values

    """
    network = target.project.network
    Ts = network.throats(target.name)
    conns = network['throat.conns']
    coords = network['pore.coords']
    return _np.mean(coords[conns], axis=1)[Ts]
