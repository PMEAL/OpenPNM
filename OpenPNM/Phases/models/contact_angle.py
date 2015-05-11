r"""
===============================================================================
Submodule -- contact_angle
===============================================================================

"""
import scipy as sp


def young(phase, sigma_sg, sigma_sl,
          surface_tension='pore.surface_tension', **kwargs):
    r"""
    Calculate contact angle using Young's equation

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
        sigma = phase[surface_tension]
        sigma = phase.interpolate_data(data=sigma)
    else:
        sigma = phase[surface_tension]
    theta = sp.arccos((sigma_sg - sigma_sl)/phase[surface_tension])
    theta = sp.rad2deg(theta)
    return theta
