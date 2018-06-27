from numpy import pi as _pi
import numpy as _np


def cylinder(target, throat_diameter='throat.diameter'):
    r"""
    Calculate throat cross-sectional area for a cylindrical throat

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_diameter : string
        Dictionary key of the throat diameter values

    """
    diams = target[throat_diameter]
    value = _pi/4*(diams)**2
    return value


def cuboid(target, throat_diameter='throat.diameter'):
    r"""
    Calculate throat cross-sectional area for a cuboid throat

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_diameter : string
        Dictionary key of the throat diameter values

    """
    diams = target[throat_diameter]
    value = (diams)**2
    return value


def equivalent_area_spherical_pores(target, throat_diameter='throat.diameter',
                                    throat_area='throat.area',
                                    pore_diameter='pore.diameter'):
    r"""
    Calculate equivalent cross-sectional area for a conduit consisting of two
    spherical pores and a constant cross-section throat. This area can be later
    used to calculate hydraulic or diffusive conductance of the conduit.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_diameter : string
        Dictionary key of the throat diameter values

    throat_area : string
        Dictionary key of the throat area values

    pore_diameter : string
        Dictionary key of the pore diameter values

    """
    network = target.project.network
    cn = network['throat.conns']
    dp = target[pore_diameter]
    d1 = dp[cn[:, 0]]
    d2 = dp[cn[:, 1]]
    dt = target[throat_diameter]
    L1 = _np.sqrt(d1**2 - dt**2) / 2            # Effective length of pore 1
    L2 = _np.sqrt(d2**2 - dt**2) / 2            # Effective length of pore 2
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    L = _np.sqrt(_np.sum((C1-C2)**2, axis=1))   # c2c distance of pores
    Lt = L - (L1+L2)                            # Effective length of throat
    return L / (2 / (d1*_pi) * _np.arctanh(2*L1/d1) + Lt / target[throat_area] +
                2 / (d2*_pi) * _np.arctanh(2*L2/d2))


def equivalent_area_circular_pores(target, throat_diameter='throat.diameter',
                                   pore_diameter='pore.diameter'):
    r"""
    Calculate equivalent cross-sectional area for a conduit consisting of two
    spherical pores and a constant cross-section throat. This area can be later
    used to calculate hydraulic or diffusive conductance of the conduit.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_diameter : string
        Dictionary key of the throat diameter values

    throat_area : string
        Dictionary key of the throat area values

    pore_diameter : string
        Dictionary key of the pore diameter values

    """
    network = target.project.network
    cn = network['throat.conns']
    dp = target[pore_diameter]
    d1 = dp[cn[:, 0]]
    d2 = dp[cn[:, 1]]
    dt = target[throat_diameter]
    L1 = _np.sqrt(d1**2 - dt**2) / 2            # Effective length of pore 1
    L2 = _np.sqrt(d2**2 - dt**2) / 2            # Effective length of pore 2
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    L = _np.sqrt(_np.sum((C1-C2)**2, axis=1))   # c2c distance of pores
    Lt = L - (L1+L2)                            # Effective length of throat
    throat_area = target[throat_diameter]
    return L / (1 / 2 * _np.arctan(2*L1/_np.sqrt(d1**2 - 4*L1**2)) +
                Lt / throat_area +
                1 / 2 * _np.arctan(2*L2/_np.sqrt(d2**2 - 4*L2**2)))
