r"""
===============================================================================
Submodule -- capillary_pressure
===============================================================================

"""

import scipy as _sp

def washburn(fluid,
             network,
             throats,
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
    fluid : OpenPNM Fluid Object
        Fluid object for the invading fluids

    Notes
    -----
    The Washburn equation is:

    .. math::
        P_c = -\frac{2\sigma(cos(\theta))}{r}

    This is the most basic approach to calculating entry pressure and is suitable for highly non-wetting invading fluids in most materials.

    """
    sigma = fluid[pore_surface_tension]
    sigma = fluid.interpolate_data(data=sigma)
    theta = fluid[pore_contact_angle]
    theta = fluid.interpolate_data(data=theta)
    r = network[throat_diameter]/2
    value = -2*sigma*_sp.cos(_sp.radians(theta))/r
    value = value[throats]
    return value

def purcell(network,
            fluid,
            throats,
            r_toroid,
            pore_surface_tension='pore.surfac_tension',
            pore_contact_angle='pore.contact_angle',
            throat_diameter='throat.diameter',
            **kwargss):
    r"""
    Computes the throat capillary entry pressure assuming the throat is a toroid.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network on which to apply the calculation
    sigma : float, array_like
        Surface tension of the invading/defending fluid pair.  Units must be consistent with the throat size values, but SI is encouraged.
    theta : float, array_like
        Contact angle formed by a droplet of the invading fluid and solid surface, measured through the defending fluid phase.  Angle must be given in degrees.
    r_toroid : float or array_like
        The radius of the solid

    Notes
    -----
    This approach accounts for the converging-diverging nature of many throat types.  Advancing the meniscus beyond the apex of the toroid requires an increase in capillary pressure beyond that for a cylindical tube of the same radius. The details of this equation are described by Mason and Morrow [1]_, and explored by Gostick [2]_ in the context of a pore network model.

    References
    ----------

    .. [1] G. Mason, N. R. Morrow, Effect of contact angle on capillary displacement curvatures in pore throats formed by spheres. J. Colloid Interface Sci. 168, 130 (1994).
    .. [2] J. Gostick, Random pore network modeling of fibrous PEMFC gas diffusion media using Voronoi and Delaunay tessellations. J. Electrochem. Soc. 160, F731 (2013).

    TODO: Triple check the accuracy of this equation
    """
    sigma = fluid[pore_surface_tension]
    sigma = fluid.interpolate_data(data=sigma)
    theta = fluid[pore_contact_angle]
    theta = fluid.interpolate_data(data=theta)
    r = network[throat_diameter]/2
    R = r_toroid
    alpha = theta - 180 + _sp.arcsin(_sp.sin(_sp.radians(theta)/(1+r/R)))
    value = (-2*sigma/r)*(_sp.cos(_sp.radians(theta - alpha))/(1 + R/r*(1-_sp.cos(_sp.radians(alpha)))))
    value = value[throats]
    return value

