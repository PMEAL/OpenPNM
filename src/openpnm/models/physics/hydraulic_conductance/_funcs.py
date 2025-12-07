import numpy as _np
from openpnm.models import _doctxt


__all__ = [
    "generic_hydraulic",
    "hagen_poiseuille",
    "valvatne_blunt",
]


@_doctxt
def generic_hydraulic(
    phase,
    pore_viscosity='pore.viscosity',
    throat_viscosity='throat.viscosity',
    size_factors='throat.hydraulic_size_factors'
):
    r"""
    Calculates the hydraulic conductance of conduits in network.

    Parameters
    ----------
    %(phase)s
    pore_viscosity : str
        %(dict_blurb)s pore viscosity
    throat_viscosity : str
        %(dict_blurb)s throat viscosity
    size_factors : str
        %(dict_blurb)s conduit hydraulic size factors

    Returns
    -------
    %(return_arr)s hydraulic conductance

    """
    network = phase.network
    conns = network.conns
    mu1, mu2 = phase[pore_viscosity][conns].T
    mut = phase[throat_viscosity]

    SF = network[size_factors]
    if isinstance(SF, dict):  # Legacy approach
        F1, Ft, F2 = SF.values()
    elif SF.ndim > 1:  # Nt-by-3 array
        F1, Ft, F2 = SF.T
    else:  # Nt array, like from network extraction predictions
        F1, Ft, F2 = _np.inf, SF, _np.inf

    g1 = F1 / mu1
    gt = Ft / mut
    g2 = F2 / mu2
    return 1 / (1/g1 + 1/gt + 1/g2)


@_doctxt
def hagen_poiseuille(
    phase,
    pore_viscosity="pore.viscosity",
    throat_viscosity="throat.viscosity",
    size_factors="throat.hydraulic_size_factors"
):
    r"""
    Calculates the hydraulic conductance of conduits in network.

    Parameters
    ----------
    %(phase)s
    pore_viscosity : str
        %(dict_blurb)s pore viscosity
    throat_viscosity : str
        %(dict_blurb)s throat viscosity
    size_factors : str
        %(dict_blurb)s conduit hydraulic size factors

    Returns
    -------
    %(return_arr)s hydraulic conductance

    """
    return generic_hydraulic(phase=phase,
                             pore_viscosity=pore_viscosity,
                             throat_viscosity=throat_viscosity,
                             size_factors=size_factors)


@_doctxt
def valvatne_blunt(
    phase,
    pore_viscosity="pore.viscosity",
    throat_viscosity="throat.viscosity",
    pore_shape_factor="pore.shape_factor",
    throat_shape_factor="throat.shape_factor",
    pore_area="pore.area",
    throat_area="throat.cross_sectional_area",
    conduit_lengths="throat.conduit_lengths",
):
    r"""
    Calculates the single phase hydraulic conductance of conduits.

    Function has been adapted for use with the Statoil imported networks
    and makes use of the shape factor in these networks to apply
    Hagen-Poiseuille flow for conduits of different shape classes:
    triangular, square and circular [2].

    Parameters
    ----------
    %(phase)s
    pore_viscosity : str
        %(dict_blurb)s pore viscosity
    throat_viscosity : str
        %(dict_blurb)s throat viscosity
    pore_shape_factor : str
        %(dict_blurb)s pore geometric shape factor
    throat_shape_factor : str
        %(dict_blurb)s throat geometric shape factor
    pore_area : str
        %(dict_blurb)s pore area
        The pore area is calculated using following formula:

        .. math::

            A_P = \frac{R_P^2}{(4 \cdot SF_P)}

        where theoratical value of pore_shape_factor in a circular tube is
        calculated using following formula:

        .. math::

            SF_P = \frac{A_P}{P_P^2} = 1/4π

    throat_area : str
        %(dict_blurb)s throat area.
        The throat area is calculated using following formula:

        .. math::

            T_A = \frac{R_T^2}{(4 \cdot SF_T)}

        where theoratical value of throat shape factor in circular tube is
        calculated using :

        .. math::

            SF_T = \frac{T_A}{T_P^2} = 1/4π

    conduit_lengths : str
        %(dict_blurb)s throat conduit lengths

    Returns
    -------
    %(return_arr)s

    References
    ----------
    [1] Valvatne, Per H., and Martin J. Blunt. "Predictive pore‐scale
    modeling of two‐phase flow in mixed wet media." Water Resources
    Research 40, no. 7 (2004).

    [2] Patzek, T. W., and D. B. Silin (2001), Shape factor and hydraulic
    conductance in noncircular capillaries I. One-phase creeping flow,
    J. Colloid Interface Sci., 236, 295–304.

    """
    network = phase.network
    conns = network["throat.conns"]
    mu_p = phase[pore_viscosity]
    mu_t = phase[throat_viscosity]

    # Fetch model parameters
    L1, Lt, L2 = network[conduit_lengths].T
    Gp = network[pore_shape_factor]
    Gt = network[throat_shape_factor]
    Ap = network[pore_area]
    At = network[throat_area]

    # Throat portions
    tri = Gt <= _np.sqrt(3) / 36.0
    circ = Gt >= 0.07
    square = ~(tri | circ)
    ntri = _np.sum(tri)
    nsquare = _np.sum(square)
    ncirc = _np.sum(circ)
    kt = _np.ones_like(Gt)
    kt[tri] = 3.0 / 5.0
    kt[square] = 0.5623
    kt[circ] = 0.5

    # Pore portions
    tri = Gp <= _np.sqrt(3) / 36.0
    circ = Gp >= 0.07
    square = ~(tri | circ)
    ntri += _np.sum(tri)
    nsquare += _np.sum(square)
    ncirc += _np.sum(circ)
    kp = _np.ones_like(Gp)
    kp[tri] = 3.0 / 5.0
    kp[square] = 0.5623
    kp[circ] = 0.5

    # Calculate conductance values
    gp = kp * Ap**2 * Gp / mu_p
    gt = kt * At**2 * Gt / mu_t
    value = L1 / gp[conns[:, 0]] + Lt / gt + L2 / gp[conns[:, 1]]

    return 1 / value
