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
    Young's equation is:
    
    .. math::
    
    \gamma_\mathrm{LG} \cos \theta_\mathrm{C} \ = \gamma_\mathrm{SG} - \gamma_\mathrm{SL}
    
    '''
    theta = sp.arccos((sigma_sg - sigma_sl)/phase[pore_surface_tension])
    theta = sp.rad2deg(theta)
    return theta