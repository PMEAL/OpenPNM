r"""
===============================================================================
Submodule -- throat_centroid
===============================================================================

"""
import numpy as _np


def pore_coords(target):
    r"""
    Calculate throat centroid values by averaging adjacent pore coordinates.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    Returns
    -------
    value : NumPy ndarray
        Array containing throat centroid values.

    """
    network = target.project.network
    Ts = network.throats(target.name)
    conns = network['throat.conns']
    coords = network['pore.coords']
    return _np.mean(coords[conns], axis=1)[Ts]
