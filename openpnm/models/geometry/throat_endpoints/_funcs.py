import numpy as _np
from openpnm.models.geometry import _geodocs


__all__ = [
    'spheres_and_cylinders',
]


@_geodocs
def spheres_and_cylinders(
    network,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter',
    throat_centroid='throat.centroid',
):
    r"""
    Computes the endpoints of throats when pores are spherical.

    The endpoints lie inside the sphere, defined by the lens formed between
    the intersection of the sphere and cylinder.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s
    %(Tcen)s

    Returns
    -------
    endpoints : ndarray
        An array containing the Nt-2-3 coordinates of the end points of each
        throat.
    """
    xyz = network.coords
    cn = network.conns
    L = network['throat.spacing']
    Dt = network[throat_diameter]
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    L1 = _np.zeros_like(L)
    L2 = _np.zeros_like(L)
    # Handle the case where Dt > Dp
    mask = Dt > D1
    L1[mask] = 0.5 * D1[mask]
    L1[~mask] = _np.sqrt(D1[~mask]**2 - Dt[~mask]**2) / 2
    mask = Dt > D2
    L2[mask] = 0.5 * D2[mask]
    L2[~mask] = _np.sqrt(D2[~mask]**2 - Dt[~mask]**2) / 2
    # Handle non-colinear pores and throat centroids
    try:
        TC = network[throat_centroid]
        LP1T = _np.linalg.norm(TC - xyz[cn[:, 0]], axis=1) + 1e-15
        LP2T = _np.linalg.norm(TC - xyz[cn[:, 1]], axis=1) + 1e-15
        unit_vec_P1T = (TC - xyz[cn[:, 0]]) / LP1T[:, None]
        unit_vec_P2T = (TC - xyz[cn[:, 1]]) / LP2T[:, None]
    except KeyError:
        unit_vec_P1T = (xyz[cn[:, 1]] - xyz[cn[:, 0]]) / L[:, None]
        unit_vec_P2T = -1 * unit_vec_P1T
    # Find throat endpoints
    EP1 = xyz[cn[:, 0]] + L1[:, None] * unit_vec_P1T
    EP2 = xyz[cn[:, 1]] + L2[:, None] * unit_vec_P2T
    # Handle throats w/ overlapping pores
    L1 = (4 * L**2 + D1**2 - D2**2) / (8 * L)
    L2 = (4 * L**2 + D2**2 - D1**2) / (8 * L)
    h = (2 * _np.sqrt(D1**2 / 4 - L1**2)).real
    overlap = L - 0.5 * (D1 + D2) < 0
    mask = overlap & (Dt < h)
    EP1[mask] = (xyz[cn[:, 0]] + L1[:, None] * unit_vec_P1T)[mask]
    EP2[mask] = (xyz[cn[:, 1]] + L2[:, None] * unit_vec_P2T)[mask]
    return {'head': EP1, 'tail': EP2}
