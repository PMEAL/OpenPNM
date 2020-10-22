import numpy as _np
from numpy import pi
from numpy import arctanh as _atanh


def spheres_and_cylinders(target,
                          pore_diameter='pore.diameter',
                          throat_diameter='throat.diameter',
                          throat_length=None,
                          conduit_lengths=None,
                          return_elements=False):
    r"""
    Compute diffusive shape coefficient for conduits of spheres and cylinders

    Parameter
    ---------
    target: OpenPNM object

    Notes
    -----
    The diffusive shape coefficient is the geometrical part of the pre-factor
    in Fick's law:

    .. math::

        n_A = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    Dt = network[throat_diameter][throats]
    # Getting areas
    A1 = (pi/4*D1**2)
    A2 = (pi/4*D2**2)
    At = (pi/4*Dt**2)
    if conduit_lengths is not None:
        L1 = network[conduit_lengths + '.pore1'][throats]
        L2 = network[conduit_lengths + '.pore2'][throats]
        Lt = network[conduit_lengths + '.throat'][throats]
    else:
        a = target[throat_diameter][throats]
        r = target[pore_diameter][cn]
        theta = _np.arcsin(_np.atleast_2d(a).T/r)
        L1, L2 = (r*_np.cos(theta)).T
        if throat_length is not None:
            Lt = target[throat_length][throats]
        else:
            C1 = network['pore.coords'][cn[:, 0]]
            C2 = network['pore.coords'][cn[:, 1]]
            L = _np.sqrt(_np.sum((C1 - C2)**2, axis=1))
            Lt = L - L1 - L2
    # Find g for half of pore 1, the throat, and half of pore 2
    g1, g2, gt = A1/L1, A2/L2, At/Lt
    # Calculate Shape factors
    # Preallocating F, SF
    # F is INTEGRAL(1/A) dx , x : 0 --> L
    F1, F2, Ft = _np.zeros((3, len(Lt)))
    SF1, SF2, SFt = _np.ones((3, len(Lt)))
    # Setting SF to 1 when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    SF1[~m1] = SF2[~m2] = SFt[~mt] = 1
    if ((_np.sum(D1 <= 2*L1) != 0) or (_np.sum(D2 <= 2*L2) != 0)):
        raise Exception('Some pores can not be modeled with ball_and_stick'
                        + 'flow shape factor. Use another model for those pores'
                        + 'with (D/L)<=2')
    # Handle the case where Dt >= Dp
    M1, M2 = [(Di <= Dt) & mi for Di, mi in zip([D1, D2], [m1, m2])]
    F1[M1] = (4*L1/(D1*Dt*pi))[M1]
    F2[M2] = (4*L2/(D2*Dt*pi))[M2]
    # Handle the rest (true balls and sticks)
    N1, N2 = [(Di > Dt) & mi for Di, mi in zip([D1, D2], [m1, m2])]
    F1[N1] = (2/(D1*pi) * _atanh(2*L1/D1))[N1]
    F2[N2] = (2/(D2*pi) * _atanh(2*L2/D2))[N2]
    Ft[mt] = (Lt/At)[mt]
    # Calculate conduit shape factors
    SF1[m1] = (L1 / (A1*F1))[m1]
    SF2[m2] = (L2 / (A2*F2))[m2]
    SFt[mt] = (Lt / (At*Ft))[mt]
    # Apply shape factors to individual g
    g1, g2, gt = g1*SF1, g2*SF2, gt*SFt
    # Ensure infinite conductance for elements with zero length
    g1[L1 == 0] = _np.inf
    g2[L2 == 0] = _np.inf
    gt[Lt == 0] = _np.inf
    if return_elements:
        vals = {'pore1': g1, 'throat': gt, 'pore2': g2}
    else:
        vals = (1/gt + 1/g1 + 1/g2)**(-1)
    return vals


def cylinders_in_series(target,
                        pore_diameter='pore.diameter',
                        throat_diameter='throat.diameter',
                        n_cylinders=5,
                        throat_length=None,
                        return_elements=False):
    r"""

    """
    network = target.network


def pyramids_and_cuboids(target,
                         pore_diameter='pore.diameter',
                         throat_diameter='throat.diameter',
                         throat_length=None,
                         return_elements=False):
    r"""

    """
    network = target.network
    R1, R2 = (target[pore_diameter][network.conns]/2).T
    Rt = target[throat_diameter]/2


def conical_frustum_and_stick(target,
                              pore_diameter='pore.diameter',
                              throat_diameter='throat.diameter',
                              pore_area='pore.area',
                              throat_area='throat.area',
                              conduit_lengths='throat.conduit_lengths',
                              return_elements=False):
    r"""
    Compute diffusive shape coefficient for conduits of spheres and cylinders
    assuming pores and throats are truncated pyramids (frustums) and constant
    cross-section cylinders (sticks).

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : string
        Dictionary key of the pore diameter values

    throat_diameter : string
        Dictionary key of the throat diameter values

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    conduit_lengths : string
        Dictionary key of the conduit lengths' values

    return_elements : Bool
        A boolean variable to whether return the values elementwise for
        pores and throat or return the final value for the conduit.

    Returns
    -------
    Either an array of diffusive_shape_coefficient or a dictionary containing the
    diffusive_shape_coefficient, which can be accessed via the dict keys 'pore1',
    'pore2', and 'throat'.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    # Get pore diameter
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    Dt = network[throat_diameter][throats]
    # Get pore/throat baseline areas (the one used in generic conductance)
    A1 = network[pore_area][cn[:, 0]]
    A2 = network[pore_area][cn[:, 1]]
    At = network[throat_area][throats]
    # Get conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    # Find g for half of pore 1, the throat, and half of pore 2
    g1, g2, gt = A1/L1, A2/L2, At/Lt
    # Preallocating F, SF
    # F is INTEGRAL(1/A) dx , x : 0 --> L
    F1, F2, Ft = _np.zeros((3, len(Lt)))
    SF1, SF2, SFt = _np.ones((3, len(Lt)))
    # Setting SF to 1 when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    SF1[~m1] = SF2[~m2] = SFt[~mt] = 1
    # Handle the rest (non-zero-length conduits)
    F1[m1] = (4*L1/(D1*Dt*pi))[m1]
    F2[m2] = (4*L2/(D2*Dt*pi))[m2]
    Ft[mt] = (Lt/At)[mt]
    # Calculate conduit shape factors
    SF1[m1] = (L1 / (A1*F1))[m1]
    SF2[m2] = (L2 / (A2*F2))[m2]
    SFt[mt] = (Lt / (At*Ft))[mt]
    # Apply shape factors to individual g
    g1, g2, gt = g1*SF1, g2*SF2, gt*SFt
    # Ensure infinite conductance for elements with zero length
    g1[L1 == 0] = _np.inf
    g2[L2 == 0] = _np.inf
    gt[Lt == 0] = _np.inf
    if return_elements:
        vals = {'pore1': g1, 'throat': gt, 'pore2': g2}
    else:
        vals = (1/gt + 1/g1 + 1/g2)**(-1)
    return vals


def cones_and_cylinders(target,
                        pore_diameter='pore.diameter',
                        throat_diameter='throat.diameter',
                        throat_length='throat.length',
                        return_elements=False):
    r"""

    """
    network = target.network
    R1, R2 = (target[pore_diameter][network.conns]/2).T
    Rt = target[throat_diameter]/2
    Lt = target[throat_length]
    L1 = Lt - R1
    L2 = Lt - R2

    alpha1 = (R1 - Rt)/L1
    beta1 = 1 / (1/(Rt**3) - 1/(R1**3))
    alpha2 = (R2-Rt)/L2
    beta2 = 1 / (1/(Rt**3) - 1/(R2**3))
    g1 = (3*alpha1*pi/8) * beta1
    g2 = (3*alpha2*pi/8) * beta2
    gt = pi*Rt**4/(8*Lt)

    if return_elements:
        g = {'pore1': g1, 'throat': gt, 'pore2': g2}
    else:
        g = (1/g1 + 1/gt + 1/g2)**-1
    return g
