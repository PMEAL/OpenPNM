
"""
module CapillaryPressure
===============================================================================

"""

import scipy as sp

def Washburn(net,fluid):
    r"""
    Computes the capillary entry pressure assuming the throat is a cylindrical tube.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network on which to apply the calculation

    fluid1 : OpenPNM Fluid Object
        Fluid object for the invading fluids

    Notes
    -----
    The Washburn equation is:

    .. math::
        P_c = -\frac{2\sigma(cos(\theta))}{r}

    This is the most basic approach to calculating entry pressure and is suitable for highly non-wetting invading fluids in most materials.

    """
    sigma = fluid.pore_conditions['surface_tension']
    sigma = fluid.interpolate_throat_conditions(net,sigma)
    theta = fluid.pore_conditions['contact_angle']
    theta = fluid.interpolate_throat_conditions(net,theta)
    vals = -4*sigma*sp.cos(sp.radians(theta))/net.throat_properties['diameter']
    fluid.throat_conditions['Pc_entry']= vals

def Purcell(net,fluid,r_toroid):
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

    """
    #This seesm to work, but I wrote it quickly and lost track of the degree-radians conversions
    """TODO:
    Triple check the accuracy of this equation
    """
    sigma = fluid.pore_conditions['surface_tension']
    sigma = fluid.interpolate_throat_conditions(net,sigma)
    theta = fluid.throat_conditions['contact_angle']
    theta = fluid.interpolate_throat_conditions(net,theta)
    r = net.throat_properties['diameter']/2
    R = r_toroid
    alpha = theta - 180 + sp.arcsin(sp.sin(sp.radians(theta)/(1+r/R)))
    vals = (-2*sigma/r)*(sp.cos(sp.radians(theta - alpha))/(1 + R/r*(1-sp.cos(sp.radians(alpha)))))
    net.throat_conditions['Pc_entry'] = vals

