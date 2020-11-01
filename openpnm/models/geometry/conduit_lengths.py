import numpy as _np
from numpy.linalg import norm as _norm

__all__ = ["spheres_and_cylinders", "circles_and_rectangles"]


def spheres_and_cylinders(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter"
):
    network = target.project.network
    throats = target.throats(target=network)
    conns = network.conns[throats]
    P1, P2 = conns.T
    X1, X2 = network.coords[conns.T]
    L_ctc = _norm(X1 - X2, axis=1)
    Dt = network[throat_diameter][throats]
    D1, D2 = network[pore_diameter][conns]
    # Handle the case where Dt > Dp
    if (Dt > D1).any() or (Dt > D2).any():
        _raise_incompatible_data()
    L1 = _np.sqrt(D1**2 - Dt**2) / 2
    L2 = _np.sqrt(D2**2 - Dt**2) / 2
    # Handle throats w/ overlapping pores
    _L1 = (4 * L_ctc**2 + D1**2 - D2**2) / (8 * L_ctc)
    mask = L_ctc - 0.5 * (D1 + D2) < 0
    L1[mask] = _L1[mask]
    L2[mask] = (L_ctc - L1)[mask]
    Lt = L_ctc - (L1 + L2)
    return _np.vstack((L1, Lt, L2)).T


def circles_and_rectangles(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter"
):
    return spheres_and_cylinders(
        target=target,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    )


# Dealing with errors and exceptions
def _raise_incompatible_data():
    raise Exception(
        "'spheres_and_cylinders' can only be applied when throat diameter is"
        " smaller than that of adjacent pores."
    )
