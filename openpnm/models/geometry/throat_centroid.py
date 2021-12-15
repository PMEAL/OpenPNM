import numpy as _np
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.dedent
def pore_coords(target):
    r"""
    Calculate throat centroid values by averaging adjacent pore coordinates

    Parameters
    ----------
    %(models.target.parameters)s

    Returns
    -------
    values : ndarray
        A numpy ndarray containing throat centroid values

    """
    network = target.project.network
    Ts = network.throats(target.name)
    conns = network['throat.conns']
    coords = network['pore.coords']
    return _np.mean(coords[conns], axis=1)[Ts]
