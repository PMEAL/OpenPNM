r"""
===============================================================================
Submodule -- contact_angle
===============================================================================

"""
import scipy as sp


def young(target, sigma_sg, sigma_sl, surface_tension='pore.surface_tension'):
    r"""
    Calculate contact angle using Young's equation

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    Notes
    -----
    Young's equation is: sigma_lg*Cos(theta) = sigma_sg - sigma_sl
    where
    sigma_lg is the liquid-gas surface tension [N/m]
    sigma_sg is the solid-gas surface tension [N/m]
    sigma_sl is the solid-liquid interfacial tension [J/m^2]
    theta is the Young contact angle [rad]

    """
    if surface_tension.split('.')[0] == 'pore':
        sigma = target[surface_tension]
        sigma = target.interpolate_data(propname=surface_tension)
    else:
        sigma = target[surface_tension]
    theta = sp.arccos((sigma_sg - sigma_sl)/target[surface_tension])
    theta = sp.rad2deg(theta)
    return theta
