r"""
===============================================================================
Submodule -- throat_centroid
===============================================================================

"""
import scipy as _sp


def pore_coords(target):
    r"""
    The average of the pore coords
    """
    network = target.project.network
    Ts = network.throats(target.name)
    conns = network['throat.conns']
    coords = network['pore.coords']
    return _sp.mean(coords[conns], axis=1)[Ts]
