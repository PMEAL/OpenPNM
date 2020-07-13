r"""

.. autofunction:: openpnm.models.physics.flow_shape_factors.ball_and_stick
.. autofunction:: openpnm.models.physics.flow_shape_factors.conical_frustum_and_stick

"""
from scipy import pi as _pi
from numpy import arctanh as _atanh
import numpy as _np
import scipy as _sp


def ball_and_stick(target, pore_area='pore.area',
                   throat_area='throat.area',
                   pore_diameter='pore.diameter',
                   throat_diameter='throat.diameter',
                   conduit_lengths='throat.conduit_lengths'):
    r"""
    Calculate conduit shape factors for hydraulic conductance, assuming
    pores and throats are spheres (balls) and constant cross-section
    cylinders (sticks).

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    pore_diameter : string
        Dictionary key of the pore diameter values

    throat_diameter : string
        Dictionary key of the throat diameter values

    conduit_lengths : string
        Dictionary key of the conduit lengths' values

    Returns
    -------
    SF : dictionary
        Dictionary containing conduit shape factors to be used in hagen-
        poiseuille hydraulic conductance model. Shape factors are accessible
        via the keys: 'pore1', 'pore2' and 'throat'.

    Notes
    -----
    (1) This model accounts for the variable cross-section area in spheres.

    (2) WARNING: This model could break if `conduit_lengths` does not
    correspond to an actual ball and stick! Example: pore length is greater
    than pore radius --> :(

    References
    ----------
    Akbari, M., Sinton, D., & Bahrami, M. (2011). Viscous flow in variable
    cross-section microchannels of arbitrary shapes. International Journal of
    Heat and Mass Transfer, 54(17-18), 3970-3978.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    # Get pore diameter
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    Dt = network[throat_diameter][throats]
    # Get conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    # Get pore/throat baseline areas (the one used in generic conductance)
    A1 = network[pore_area][cn[:, 0]]
    A2 = network[pore_area][cn[:, 1]]
    At = network[throat_area][throats]
    # Preallocating F, SF
    # F is INTEGRAL(1/A^2) dx , x : 0 --> L
    F1, F2, Ft = _np.zeros((3, len(Lt)))
    SF1, SF2, SFt = _np.ones((3, len(Lt)))
    # Setting SF to 1 when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    SF1[~m1] = SF2[~m2] = SFt[~mt] = 1
    # Handle the case where Dt >= Dp
    M1, M2 = [(Di <= Dt) & mi for Di, mi in zip([D1, D2], [m1, m2])]
    F1[M1] = 16/3 * (L1*(D1**2 + D1*Dt + Dt**2) / (D1**3 * Dt**3 * _pi**2))[M1]
    F2[M2] = 16/3 * (L2*(D2**2 + D2*Dt + Dt**2) / (D2**3 * Dt**3 * _pi**2))[M2]
    # Handle the rest (true balls and sticks)
    N1, N2 = [(Di > Dt) & mi for Di, mi in zip([D1, D2], [m1, m2])]
    F1[N1] = (4/(D1**3*_pi**2) * ((2*D1*L1) / (D1**2-4*L1**2) + _atanh(2*L1/D1)))[N1]
    F2[N2] = (4/(D2**3*_pi**2) * ((2*D2*L2) / (D2**2-4*L2**2) + _atanh(2*L2/D2)))[N2]
    Ft[mt] = (Lt / At**2)[mt]
    # Calculate conduit shape factors
    SF1[m1] = (L1 / (A1**2 * F1))[m1]
    SF2[m2] = (L2 / (A2**2 * F2))[m2]
    SFt[mt] = (Lt / (At**2 * Ft))[mt]
    return {'pore1': SF1, 'throat': SFt, 'pore2': SF2}


def conical_frustum_and_stick(target, pore_area='pore.area',
                              throat_area='throat.area',
                              pore_diameter='pore.diameter',
                              throat_diameter='throat.diameter',
                              conduit_lengths='throat.conduit_lengths'):
    r"""
    Calculate conduit shape factors for hydraulic conductance, assuming
    pores and throats are truncated pyramids (frustums) and constant
    cross-section cylinders (sticks).

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    pore_diameter : string
        Dictionary key of the pore diameter values

    throat_diameter : string
        Dictionary key of the throat diameter values

    conduit_lengths : string
        Dictionary key of the conduit lengths' values

    Returns
    -------
    SF : dictionary
        Dictionary containing conduit shape factors to be used in hagen-
        poiseuille hydraulic conductance model. Shape factors are accessible
        via the keys: 'pore1', 'pore2' and 'throat'.

    References
    ----------
    Akbari, M., Sinton, D., & Bahrami, M. (2011). Viscous flow in variable
    cross-section microchannels of arbitrary shapes. International Journal of
    Heat and Mass Transfer, 54(17-18), 3970-3978.

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
    # Preallocating F, SF
    # F is INTEGRAL(1/A^2) dx , x : 0 --> L
    F1, F2, Ft = _np.zeros((3, len(Lt)))
    SF1, SF2, SFt = _np.ones((3, len(Lt)))
    # Setting SF to 1 when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    SF1[~m1] = SF2[~m2] = SFt[~mt] = 1
    # Calculate integral of 1/A^2
    F1[m1] = 16/3 * (L1*(D1**2 + D1*Dt + Dt**2) / (D1**3 * Dt**3 * _pi**2))[m1]
    F2[m2] = 16/3 * (L2*(D2**2 + D2*Dt + Dt**2) / (D2**3 * Dt**3 * _pi**2))[m2]
    Ft[mt] = (Lt/At**2)[mt]
    # Calculate conduit shape factors
    SF1[m1] = (L1 / (A1**2 * F1))[m1]
    SF2[m2] = (L2 / (A2**2 * F2))[m2]
    SFt[mt] = (Lt / (At**2 * Ft))[mt]
    return {'pore1': SF1, 'throat': SFt, 'pore2': SF2}


def ball_and_stick_2D(target, pore_area='pore.area',
                      throat_area='throat.area',
                      pore_diameter='pore.diameter',
                      throat_diameter='throat.diameter',
                      conduit_lengths='throat.conduit_lengths'):
    r"""
    Calculate conduit shape factors for hydraulic conductance, assuming
    pores and throats are circles (balls) and rectangles (sticks).

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    pore_diameter : string
        Dictionary key of the pore diameter values

    throat_diameter : string
        Dictionary key of the throat diameter values

    conduit_lengths : string
        Dictionary key of the conduit lengths' values

    Returns
    -------
    SF : dictionary
        Dictionary containing conduit shape factors to be used in hagen-
        poiseuille hydraulic conductance model. Shape factors are accessible
        via the keys: 'pore1', 'pore2' and 'throat'.

    Notes
    -----
    (1) This model accounts for the variable cross-section area in spheres.

    (2) WARNING: This model could break if `conduit_lengths` does not
    correspond to an actual ball and stick! Example: pore length is greater
    than pore radius --> :(

    References
    ----------
    Akbari, M., Sinton, D., & Bahrami, M. (2011). Viscous flow in variable
    cross-section microchannels of arbitrary shapes. International Journal of
    Heat and Mass Transfer, 54(17-18), 3970-3978.

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
    A1 = network[pore_area][cn[:, 0]]
    A2 = network[pore_area][cn[:, 1]]
    At = network[throat_area][throats]
    # Preallocating F, SF
    # F is INTEGRAL(1/A^2) dx , x : 0 --> L
    F1, F2, Ft = _np.zeros((3, len(Lt)))
    SF1, SF2, SFt = _np.ones((3, len(Lt)))
    # Setting SF to 1 when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    SF1[~m1] = SF2[~m2] = SFt[~mt] = 1
    F1[m1] = (_atanh(2*L1/D1) / (2*D1))[m1]
    F2[m2] = (_atanh(2*L2/D2) / (2*D2))[m2]
    Ft[mt] = (Lt / At**2)[mt]
    # Calculate conduit shape factors
    SF1[m1] = (L1 / (A1**2 * F1))[m1]
    SF2[m2] = (L2 / (A2**2 * F2))[m2]
    SFt[mt] = (Lt / (At**2 * Ft))[mt]
    return {'pore1': SF1, 'throat': SFt, 'pore2': SF2}
