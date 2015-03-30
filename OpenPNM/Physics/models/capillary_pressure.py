r"""
===============================================================================
Submodule -- capillary_pressure
===============================================================================

"""

import scipy as _sp

def washburn(physics,
             phase,
             network,
             pore_surface_tension='pore.surface_tension',
             pore_contact_angle='pore.contact_angle',
             throat_diameter='throat.diameter',
             **kwargs):
    r"""
    Computes the capillary entry pressure assuming the throat is a cylindrical tube.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network on which to apply the calculation
    phase : OpenPNM Phase Object
        Phase object for the invading phases

    Notes
    -----
    The Washburn equation is:

    .. math::
        P_c = -\frac{2\sigma(cos(\theta))}{r}

    This is the most basic approach to calculating entry pressure and is 
    suitable for highly non-wetting invading phases in most materials.

    """
    sigma = phase[pore_surface_tension]
    sigma = phase.interpolate_data(data=sigma)
    theta = phase[pore_contact_angle]
    theta = phase.interpolate_data(data=theta)
    r = network[throat_diameter]/2
    value = -2*sigma*_sp.cos(_sp.radians(theta))/r
    value = value[phase.throats(physics.name)]
    return value

def purcell(physics,
            phase,
            network,
            r_toroid,
            pore_surface_tension='pore.surface_tension',
            pore_contact_angle='pore.contact_angle',
            throat_diameter='throat.diameter',
            **kwargs):
    r"""
    Computes the throat capillary entry pressure assuming the throat is a toroid.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network on which to apply the calculation
    sigma : float, array_like
        Surface tension of the invading/defending phase pair.  Units must be 
        consistent with the throat size values, but SI is encouraged.
    theta : float, array_like
        Contact angle formed by a droplet of the invading phase and solid 
        surface, measured through the defending phase phase.  Angle must be given in degrees.
    r_toroid : float or array_like
        The radius of the solid

    Notes
    -----
    This approach accounts for the converging-diverging nature of many throat 
    types.  Advancing the meniscus beyond the apex of the toroid requires an 
    increase in capillary pressure beyond that for a cylindical tube of the 
    same radius. The details of this equation are described by Mason and 
    Morrow [1]_, and explored by Gostick [2]_ in the context of a pore network 
    model.

    References
    ----------

    .. [1] G. Mason, N. R. Morrow, Effect of contact angle on capillary displacement curvatures in pore throats formed by spheres. J. Colloid Interface Sci. 168, 130 (1994).
    .. [2] J. Gostick, Random pore network modeling of fibrous PEMFC gas diffusion media using Voronoi and Delaunay tessellations. J. Electrochem. Soc. 160, F731 (2013).

    TODO: Triple check the accuracy of this equation
    """

    sigma = phase[pore_surface_tension]
    sigma = phase.interpolate_data(data=sigma)
    theta = phase[pore_contact_angle]
    theta = phase.interpolate_data(data=theta)
    r = network[throat_diameter]/2
    R = r_toroid
    alpha = theta - 180 + _sp.arcsin(_sp.sin(_sp.radians(theta)/(1+r/R)))
    value = (-2*sigma/r)*(_sp.cos(_sp.radians(theta - alpha))/(1 + R/r*(1-_sp.cos(_sp.radians(alpha)))))
    value = value[phase.throats(physics.name)]
    return value
    
def purcell_pore(physics,
            phase,
            network,
            r_toroid,
            pore_surface_tension='pore.surface_tension',
            pore_contact_angle='pore.contact_angle',
            pore_diameter='pore.diameter',
            **kwargs):
    r"""
    Computes the pore capillary entry pressure assuming the pore is a toroid.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network on which to apply the calculation
    sigma : float, array_like
        Surface tension of the invading/defending phase pair.  Units must be 
        consistent with the throat size values, but SI is encouraged.
    theta : float, array_like
        Contact angle formed by a droplet of the invading phase and solid 
        surface, measured through the defending phase phase.  Angle must be given in degrees.
    r_toroid : float or array_like
        The radius of the solid

    Notes
    -----
    This approach accounts for the converging-diverging nature of many throat 
    types.  Advancing the meniscus beyond the apex of the toroid requires an 
    increase in capillary pressure beyond that for a cylindical tube of the 
    same radius. The details of this equation are described by Mason and 
    Morrow [1]_, and explored by Gostick [2]_ in the context of a pore network 
    model.

    References
    ----------

    .. [1] G. Mason, N. R. Morrow, Effect of contact angle on capillary displacement curvatures in pore throats formed by spheres. J. Colloid Interface Sci. 168, 130 (1994).
    .. [2] J. Gostick, Random pore network modeling of fibrous PEMFC gas diffusion media using Voronoi and Delaunay tessellations. J. Electrochem. Soc. 160, F731 (2013).

    TODO: Triple check the accuracy of this equation
    """

    sigma = phase[pore_surface_tension]
    theta = phase[pore_contact_angle]
    r = network[pore_diameter]/2
    R = r_toroid
    alpha = theta - 180 + _sp.arcsin(_sp.sin(_sp.radians(theta)/(1+r/R)))
    value = (-2*sigma/r)*(_sp.cos(_sp.radians(theta - alpha))/(1 + R/r*(1-_sp.cos(_sp.radians(alpha)))))
    value = value[phase.pores(physics.name)]
    return value