import openpnm.models as mods
import numpy as _np
from numpy import pi as _pi
from numpy import arctanh as _atanh
from numpy import sqrt as _sqrt


def spheres_and_cylinders(target,
                          pore_diameter='pore.diameter',
                          throat_diameter='throat.diameter',
                          throat_length=None,
                          conduit_lengths=None):
    r"""
    Compute diffusive shape coefficient for conduits of spheres and cylinders

    Parameter
    ---------
    target: OpenPNM Geometry object


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
    L1 = network[conduit_lengths + '.pore1'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]

    # F is Integral(1/A) dx , x : 0 --> L
    if ((_np.sum(D1 <= 2*L1) != 0) or (_np.sum(D2 <= 2*L2) != 0)):
        raise Exception('Some throats are too short, add spherical_pores endpoint model')
    F1 = 2/(D1*_pi) * _atanh(2*L1/D1)
    F2 = 2/(D2*_pi) * _atanh(2*L2/D2)
    Ft = Lt / (_pi/4*Dt**2)
    g1, g2, gt = 1/F1, 1/F2, 1/Ft
    vals = {'pore1': g1, 'throat': gt, 'pore2': g2}
    return vals


def spheres_and_cylinders_2D(target,
                             pore_diameter='pore.diameter',
                             throat_diameter='throat.diameter',
                             throat_length=None,
                             conduit_lengths=None,
                             return_elements=False):
    r"""
    Compute diffusive shape coefficient for conduits of spheres and cylinders
    in 2D

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
    # Get pore diameter
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    # Get conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    # Get pore/throat baseline areas (the one used in generic conductance)
    A1 = D1
    A2 = D2
    At = network[throat_diameter][throats]
    # Find g for half of pore 1, the throat, and half of pore 2
    g1, g2, gt = A1/L1, A2/L2, At/Lt
    # Apply shape factors to individual g
    SF = _SF_spheres_and_cylinders_2D(L1, L2, Lt, D1, D2, A1, A2, At)
    SF1, SF2, SFt = SF['pore1'], SF['pore2'], SF['throat']
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
                        pore_area='pore.area',
                        throat_area='throat.area',
                        conduit_lengths='throat.conduit_lengths',
                        return_elements=False):
    r"""
    Compute diffusive shape coefficient assuming pores are truncated pyramids
    and throats are cylinders.

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
    # Apply shape factors to individual g
    mod = mods.geometry.diffusive_shape_factors.conical_frustum_and_stick
    SF = mod(target=target, pore_area=pore_area,
             throat_area=throat_area,
             pore_diameter=pore_diameter,
             throat_diameter=throat_diameter,
             conduit_lengths=conduit_lengths)
    # SF = _SF_cones_and_cylinders(L1, L2, Lt, D1, D2, Dt, A1, A2, At)
    SF1, SF2, SFt = SF['pore1'], SF['pore2'], SF['throat']
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


def pyramids_and_cuboids(target,
                         pore_diameter='pore.diameter',
                         throat_diameter='throat.diameter',
                         return_elements=False,
                         ):
    r"""

    """
    pass


def cubes_and_cuboids(target,
                      pore_diameter='pore.diameter',
                      throat_diameter='throat.diameter',
                      pore_aspect=[1, 1, 1],
                      throat_aspect=[1, 1, 1],
                      return_elements=False
                      ):
    r"""

    """
    pass


def intersecting_cones(target,
                       pore_diameter='pore.diameter',
                       throat_diameter='throat.diameter',
                       midpoint=None,
                       return_elements=False
                       ):
    r"""

    """
    network = target.network
    R1, R2 = (target[pore_diameter][network.conns]/2).T
    Rt = target[throat_diameter]/2
    Lt = 0
    L1 = Lt - R1
    L2 = Lt - R2

    alpha1 = (R1 - Rt)/L1
    beta1 = 1 / (1/(Rt**3) - 1/(R1**3))
    alpha2 = (R2-Rt)/L2
    beta2 = 1 / (1/(Rt**3) - 1/(R2**3))
    g1 = (3*alpha1*_pi/8) * beta1
    g2 = (3*alpha2*_pi/8) * beta2
    gt = _pi*Rt**4/(8*Lt)

    if return_elements:
        g = {'pore1': g1, 'throat': gt, 'pore2': g2}
    else:
        g = (1/g1 + 1/gt + 1/g2)**-1
    return g


def intersecting_pyramids(target,
                          pore_diameter='pore.diameter',
                          throat_diameter='throat.diameter',
                          midpoint=None,
                          return_elements=False):
    r"""
    """
    pass


def ncylinders_in_series(target,
                         pore_diameter='pore.equivalent_diameter',
                         throat_diameter='throat.equivalent_diameter',
                         throat_length='throat.length',
                         n=5, return_elements=True):
    r"""
    Computes the shape coefficient of pores as N cylinders in series

    Parameters
    ----------
    target : OpenPNM Geometry object
        The object to which this model applies
    pore_diameter : str
        Dictionary key pointing to the pore diameter values
    throat_diameter : str
        Dictionary key pointing to the throat diameter values
    throat_length : str
    """
    project = target.project
    network = project.network
    P12 = network['throat.conns']
    Pdia1, Pdia2 = network[pore_diameter][P12].T
    Tdia = network[throat_diameter]
    # Ensure throats are never bigger than connected pores
    Tdia = _np.minimum(Tdia, 0.99*_np.minimum(Pdia1, Pdia2))
    Plen1 = Pdia1/2*(_np.cos(_np.arcsin(Tdia/Pdia1)))
    Plen2 = Pdia2/2*(_np.cos(_np.arcsin(Tdia/Pdia2)))
    Lcyl1 = _np.linspace(0, Plen1, num=n, endpoint=False)
    Lcyl2 = _np.linspace(0, Plen2, num=n, endpoint=False)
    Rcyl1 = Pdia1/2*_np.sin(_np.arccos(Lcyl1/(Pdia1/2)))
    Rcyl2 = Pdia2/2*_np.sin(_np.arccos(Lcyl2/(Pdia2/2)))
    gtemp = (_pi*Rcyl1**2/(Plen1/n)).T
    g1 = 1/_np.sum(1/gtemp, axis=1)
    gtemp = (_pi*Rcyl2**2/(Plen2/n)).T
    g2 = 1/_np.sum(1/gtemp, axis=1)
    Tlen = network[throat_length]
    gt = (_pi*(Tdia/2)**2/(Tlen)).T
    if return_elements:
        result = {'pore1': g1, 'throat': gt, 'pore2': g2}
    else:
        result = 1/(1/g1 + 1/gt + 1/g2)
    return result


def _SF_spheres_and_cylinders(L1, L2, Lt, D1, D2, Dt, A1, A2, At):
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
    F1[M1] = (4*L1/(D1*Dt*_pi))[M1]
    F2[M2] = (4*L2/(D2*Dt*_pi))[M2]
    # Handle the rest (true balls and sticks)
    N1, N2 = [(Di > Dt) & mi for Di, mi in zip([D1, D2], [m1, m2])]
    F1[N1] = (2/(D1*_pi) * _atanh(2*L1/D1))[N1]
    F2[N2] = (2/(D2*_pi) * _atanh(2*L2/D2))[N2]
    Ft[mt] = (Lt/At)[mt]
    # Calculate conduit shape factors
    SF1[m1] = (L1 / (A1*F1))[m1]
    SF2[m2] = (L2 / (A2*F2))[m2]
    SFt[mt] = (Lt / (At*Ft))[mt]
    return {'pore1': SF1, 'throat': SFt, 'pore2': SF2}


def _SF_spheres_and_cylinders_2D(L1, L2, Lt, D1, D2, A1, A2, At):
    # Preallocating F, SF
    # F is INTEGRAL(1/A) dx , x : 0 --> L
    F1, F2, Ft = _np.zeros((3, len(Lt)))
    SF1, SF2, SFt = _np.ones((3, len(Lt)))
    # Setting SF to 1 when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    SF1[~m1] = SF2[~m2] = SFt[~mt] = 1
    F1[m1] = (0.5 * _atanh(2*L1/_sqrt(D1**2 - 4*L1**2)))[m1]
    F2[m2] = (0.5 * _atanh(2*L2/_sqrt(D2**2 - 4*L2**2)))[m2]
    Ft[mt] = (Lt/At)[mt]
    # Calculate conduit shape factors
    SF1[m1] = (L1 / (A1*F1))[m1]
    SF2[m2] = (L2 / (A2*F2))[m2]
    SFt[mt] = (Lt / (At*Ft))[mt]
    return {'pore1': SF1, 'throat': SFt, 'pore2': SF2}


def _SF_cones_and_cylinders(L1, L2, Lt, D1, D2, Dt, A1, A2, At):
    # Preallocating F, SF
    # F is INTEGRAL(1/A) dx , x : 0 --> L
    F1, F2, Ft = _np.zeros((3, len(Lt)))
    SF1, SF2, SFt = _np.ones((3, len(Lt)))
    # Setting SF to 1 when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    SF1[~m1] = SF2[~m2] = SFt[~mt] = 1
    # Handle the rest (non-zero-length conduits)
    F1[m1] = (4*L1/(D1*Dt*_pi))[m1]
    F2[m2] = (4*L2/(D2*Dt*_pi))[m2]
    Ft[mt] = (Lt/At)[mt]
    # Calculate conduit shape factors
    SF1[m1] = (L1 / (A1*F1))[m1]
    SF2[m2] = (L2 / (A2*F2))[m2]
    SFt[mt] = (Lt / (At*Ft))[mt]
    return {'pore1': SF1, 'throat': SFt, 'pore2': SF2}
