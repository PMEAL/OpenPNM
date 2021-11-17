r"""
Pore-scale models for calculating hydraulic conductance of conduits.
"""
import numpy as _np

__all__ = [
    "generic_hydraulic",
    "hagen_poiseuille",
    "hagen_poiseuille_power_law",
    "valvatne_blunt",
]


def generic_hydraulic(target,
                      pore_viscosity='pore.viscosity',
                      throat_viscosity='throat.viscosity',
                      size_factors='throat.hydraulic_size_factors'):
    r"""
    Calculates the hydraulic conductance of conduits in network.

    A conduit is defined as (1/2 pore - full throat - 1/2 pore).

    Parameters
    ----------
    target : _GenericPhysics
        Physics object with which this model is associated.
    pore_viscosity : str
        Dictionary key of the pore viscosity values.
    throat_viscosity : str
        Dictionary key of the throat viscosity values.
    size_factors: str
        Dictionary key of the conduit hydraulic size factors' values.

    Returns
    -------
    ndarray
        Array containing hydraulic conductance values for conduits in the
        geometry attached to the given physics object.

    """
    network = target.network
    throats = target.throats(to_global=True)
    conns = network.conns[throats]
    phase = target.project.find_phase(target)
    F = network[size_factors]
    mu1, mu2 = phase[pore_viscosity][conns].T
    mut = phase[throat_viscosity][throats]

    if isinstance(F, dict):
        g1 = F[f"{size_factors}.pore1"][throats] / mu1
        gt = F[f"{size_factors}.throat"][throats] / mut
        g2 = F[f"{size_factors}.pore2"][throats] / mu2
        return 1 / (1/g1 + 1/gt + 1/g2)
    return F[throats] / mut


def hagen_poiseuille(
    target,
    pore_viscosity="pore.viscosity",
    throat_viscosity="throat.viscosity",
    size_factors="throat.hydraulic_size_factors"
):
    r"""
    Calculates the hydraulic conductance of conduits in network.

    A conduit is defined as (1/2 pore - full throat - 1/2 pore).

    Parameters
    ----------
    target : _GenericPhysics
        Physics object with which this model is associated.
    pore_viscosity : str
        Dictionary key of the pore viscosity values.
    throat_viscosity : str
        Dictionary key of the throat viscosity values.
    size_factors: str
        Dictionary key of the conduit size factors' values.

    Returns
    -------
    g : ndarray
        Array containing hydraulic conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This function requires that all the necessary phase properties already
    be calculated.

    """
    return generic_hydraulic(target=target,
                             pore_viscosity=pore_viscosity,
                             throat_viscosity=throat_viscosity,
                             size_factors=size_factors)


def hagen_poiseuille_power_law(
    target,
    pore_area="pore.area",
    throat_area="throat.area",
    pore_viscosity_min="pore.viscosity_min",
    throat_viscosity_min="throat.viscosity_min",
    pore_viscosity_max="pore.viscosity_max",
    throat_viscosity_max="throat.viscosity_max",
    conduit_lengths="throat.conduit_lengths",
    size_factors="throat.hydraulic_size_factors",
    pore_consistency="pore.consistency",
    throat_consistency="throat.consistency",
    pore_flow_index="pore.flow_index",
    throat_flow_index="throat.flow_index",
    pore_pressure="pore.pressure",
):
    r"""
    Calculates the hydraulic conductance of conduits in network, assuming
    a non Newtonian fluid whose viscosity obeys a power law.

    A conduit is defined as (1/2 pore - full throat - 1/2 pore).

    Parameters
    ----------
    target : _GenericPhysics
        Physics object with which this model is associated.
    pore_area : str
        Dictionary key of the pore area values.
    throat_area : str
        Dictionary key of the throat area values.
    pore_viscosity_min : str
        Dictionary key of the pore minimum viscosity values.
    throat_viscosity_min : str
        Dictionary key of the throat minimum viscosity values.
    pore_viscosity_max : str
        Dictionary key of the pore maximum viscosity values.
    throat_viscosity_max : str
        Dictionary key of the throat maximum viscosity values.
    conduit_lengths : str
        Dictionary key of the conduit lengths' values.
    size_factors: str
        Dictionary key of the conduit size factors' values.
    pore_consistency : str
        Dictionary key of the pore fluid consistency values.
    throat_consistency : str
        Dictionary key of the throat fluid consistency values.
    pore_flow_index : str
        Dictionary key of the pore fluid flow index values.
    throat_flow_index : str
        Dictionary key of the throat fluid flow index values.
    pore_pressure : str
        Dictionary key of the pore pressure values.

    Returns
    -------
    g : ndarray
        Array containing hydraulic conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This function requires that all the necessary phase properties already
    be calculated.

    """
    # Fetch openpnm objects
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network.conns[throats]

    # Fetch model parameters
    A1, A2 = network[pore_area][cn].T
    At = network[throat_area][throats]
    L1 = network[conduit_lengths + ".pore1"][throats]
    L2 = network[conduit_lengths + ".pore2"][throats]
    Lt = network[conduit_lengths + ".throat"][throats]
    P = phase[pore_pressure]
    Pt = phase.interpolate_data(propname=pore_pressure)[throats]

    mu_min = phase[pore_viscosity_min][cn]
    mu_mint = phase[throat_viscosity_min][throats]
    mu_max = phase[pore_viscosity_max][cn]
    mu_maxt = phase[throat_viscosity_max][throats]

    C1, C2 = phase[pore_consistency][cn].T
    Ct = phase[throat_consistency][throats]
    n1, n2 = phase[pore_flow_index][cn].T
    nt = phase[throat_flow_index][throats]

    # Pressure differences dP
    dP1 = _np.absolute(P[cn[:, 0]] - Pt).clip(min=1e-20)
    dP2 = _np.absolute(P[cn[:, 1]] - Pt).clip(min=1e-20)
    dPt = _np.absolute(_np.diff(P[cn], axis=1).squeeze()).clip(min=1e-20)

    # Calculate apparent viscosities
    mu1 = (dP1**(1-1/n1) * C1**(1/n1)) / ((4*n1 / (3*n1+1)) * (2*L1 / ((A1/_np.pi)**0.5))**(1-1/n1))
    mu2 = (dP2**(1-1/n2) * C2**(1/n2)) / ((4*n2 / (3*n2+1)) * (2*L2 / ((A2/_np.pi)**0.5))**(1-1/n2))
    mut = (dPt**(1-1/nt) * Ct**(1/nt)) / ((4*nt / (3*nt+1)) * (2*Lt / ((At/_np.pi)**0.5))**(1-1/nt))

    # Bound the apparent viscosities
    mu1 = _np.minimum(_np.maximum(mu1, mu_min[:, 0]), mu_max[:, 0])
    mu2 = _np.minimum(_np.maximum(mu2, mu_min[:, 1]), mu_max[:, 1])
    mut = _np.minimum(_np.maximum(mut, mu_mint), mu_maxt)

    F = network[size_factors]
    if isinstance(F, dict):
        g1 = F[f"{size_factors}.pore1"][throats] / mu1
        gt = F[f"{size_factors}.throat"][throats] / mut
        g2 = F[f"{size_factors}.pore2"][throats] / mu2
        return 1 / (1 / g1 + 1 / gt + 1 / g2)
    return F[throats] / mut


def valvatne_blunt(
    target,
    pore_viscosity="pore.viscosity",
    throat_viscosity="throat.viscosity",
    pore_shape_factor="pore.shape_factor",
    throat_shape_factor="throat.shape_factor",
    pore_area="pore.area",
    throat_area="throat.area",
    conduit_lengths="throat.conduit_lengths",
):
    r"""
    Calculate the single phase hydraulic conductance of conduits in network,
    where a conduit is (1/2 pore - full throat - 1/2 pore) according to [1].
    Function has been adapted for use with the Statoil imported networks and
    makes use of the shape factor in these networks to apply Hagen-Poiseuille
    flow for conduits of different shape classes: Triangular, Square and
    Circular [2].

    Parameters
    ----------
    target : GenericPhysics
        Physics object with which this model is associated.
    pore_viscosity : str
        Dictionary key of the pore viscosity values.
    throat_viscosity : str
        Dictionary key of the throat viscosity values.
    pore_shape_factor : str
        Dictionary key of the pore geometric shape factor values.
    throat_shape_factor : str
        Dictionary key of the throat geometric shape factor values.
    pore_area : str
        Dictionary key of the pore area values. The pore area is
        calculated using following formula:
            pore_area = (pore_radius ** 2) / (4 * pore_shape_factor)
        Where theoratical value of pore_shape_factor in circular tube is
        calculated using following formula:
            pore_shape_factor = pore_area / perimeter **2 = 1/4π
    throat_area : str
        Dictionary key of the throat area values. The throat area is
        calculated using following formula:
            throat_area = (throat_radius ** 2) / (4 * throat_shape_factor)
        Where theoratical value of throat_shape_factor in circular tube is
        calculated using following formula:
            throat_shape_factor = throat_area / perimeter ** 2 = 1/4π
    conduit_lengths : str
        Dictionary key of the throat conduit lengths.

    Returns
    -------
    g : ndarray
        Array containing hydraulic conductance values for conduits in the
        geometry attached to the given physics object.

    References
    ----------
    [1] Valvatne, Per H., and Martin J. Blunt. "Predictive pore‐scale
    modeling of two‐phase flow in mixed wet media." Water Resources
    Research 40, no. 7 (2004).

    [2] Patzek, T. W., and D. B. Silin (2001), Shape factor and hydraulic
    conductance in noncircular capillaries I. One-phase creeping flow,
    J. Colloid Interface Sci., 236, 295–304.

    """
    # Fetch openpnm objects
    network = target.network
    conns = network["throat.conns"]
    mu_p = target[pore_viscosity]
    mu_t = target[throat_viscosity]

    # Fetch model parameters
    L1 = network[conduit_lengths + ".pore1"]
    Lt = network[conduit_lengths + ".throat"]
    L2 = network[conduit_lengths + ".pore2"]
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
