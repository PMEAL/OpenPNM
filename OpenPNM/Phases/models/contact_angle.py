r"""
===============================================================================
Submodule -- contact_angle
===============================================================================

"""
import scipy as sp

def young(phase,
          sigma_sg,
          sigma_sl,
          pore_surface_tension='pore.surface_tension',
          **kwargs):
    r'''
    Calculate contact angle using Young's equation
    
    Notes
    -----
    Young's equation is: sigma_lg*Cos(theta)=sigma_sg - sigma_sl
    where
    sigma_lg is the liquid-gas surface tension [N/m]
    sigma_sg is the solid-gas surface tension [N/m]
    sigma_sl is the solid-liquid interfacial tension [J/m^2]
    theta is the Young contact angle [rad]
           
    '''
    theta = sp.arccos((sigma_sg - sigma_sl)/phase[pore_surface_tension])
    theta = sp.rad2deg(theta)
    return theta