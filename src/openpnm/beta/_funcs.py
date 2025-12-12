import numpy as np


__all__ = [
    'relax_points',
]


def relax_points(network, iterations=3, f=0.5):
    r"""
    Performs a relaxation of Delaunay points based on distance weighted average
    of connected Voronoi points

    Parameters
    ----------
    network : dict
        OpenPNM network dictionary
    iterations : int
        The number of time to iteratively apply the relaxation
    f : scalar
        A damping factor to limit the how much the point actually moves towards
        new computed point, to prevent overshooting.

    Returns
    -------
    coords : ndarray
        The Np-by-3 array of pore coordinates after relaxation.
    """
    if isinstance(f, int):
        f = [f]*iterations
    Ts = network.throats('interconnect')
    Ps = network.pores('delaunay')
    conns = network.conns[Ts]
    mask1 = network['pore.delaunay'][conns[:, 0]]
    mask2 = network['pore.delaunay'][conns[:, 1]]
    coords = np.copy(network.coords)
    for i in range(iterations):
        new_coords = np.zeros((len(Ps), 3), dtype=float)
        d = np.sqrt(np.sum((coords[conns[:, 0]] - coords[conns[:, 1]])**2, axis=1))
        weights = np.zeros(len(Ps), dtype=float)
        for ax in range(coords.shape[1]):
            np.add.at(new_coords[:, ax],
                      conns[:, 0][mask1],
                      d[mask1]*coords[conns[:, 1][mask1], ax])
            np.add.at(new_coords[:, ax],
                      conns[:, 1][mask2],
                      d[mask2]*coords[conns[:, 0][mask2], ax])
        np.add.at(weights, conns[:, 0][mask1], d[mask1])
        np.add.at(weights, conns[:, 1][mask2], d[mask2])
        delta = ((new_coords.T/weights).T - coords[Ps, :])*f[i]
        coords[Ps, :] += delta
    return coords
