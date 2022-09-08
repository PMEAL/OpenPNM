import numpy as _np
from openpnm.models.geometry import _geodocs


__all__ = ["pore_coords"]


@_geodocs
def pore_coords(
    network
):
    r"""
    Calculate throat centroid values by averaging adjacent pore coordinates

    Parameters
    ----------
    %(network)s

    Returns
    -------
    values : ndarray
        A numpy ndarray containing throat centroid values

    """
    conns = network['throat.conns']
    coords = network['pore.coords']
    return _np.mean(coords[conns], axis=1)
