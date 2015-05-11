r"""
===============================================================================
Submodule -- capillary_pressure
===============================================================================

"""

import scipy as _sp


def washburn(physics, phase, network, surface_tension='pore.surface_tension',
             contact_angle='pore.contact_angle', throat_diameter='throat.diameter',
             **kwargs):
    r"""
    Computes the capillary entry pressure assuming the throat is a cylindrical tube.

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network object is
    phase : OpenPNM Phase Object
        Phase object for the invading phases containing the surface tension and
        contact angle values.
    sigma : dict key (string)
        The dictionary key containing the surface tension values to be used. If
        a pore property is given, it is interpolated to a throat list.
    theta : dict key (string)
        The dictionary key containing the contact angle values to be used. If
        a pore property is given, it is interpolated to a throat list.
    throat_diameter : dict key (string)
        The dictionary key containing the throat diameter values to be used.

    Notes
    -----
    The Washburn equation is:

    .. math::
        P_c = -\frac{2\sigma(cos(\theta))}{r}

    This is the most basic approach to calculating entry pressure and is
    suitable for highly non-wetting invading phases in most materials.

    """
    if surface_tension.split('.')[0] == 'pore':
        sigma = phase[surface_tension]
        sigma = phase.interpolate_data(data=sigma)
    else:
        sigma = phase[surface_tension]
    if contact_angle.split('.')[0] == 'pore':
        theta = phase[contact_angle]
        theta = phase.interpolate_data(data=theta)
    else:
        theta = phase[contact_angle]
    r = network[throat_diameter]/2
    value = -2*sigma*_sp.cos(_sp.radians(theta))/r
    value = value[phase.throats(physics.name)]
    return value


def purcell(physics, phase, network, r_toroid,
            surface_tension='pore.surface_tension',
            contact_angle='pore.contact_angle',
            throat_diameter='throat.diameter', **kwargs):
    r"""
    Computes the throat capillary entry pressure assuming the throat is a toroid.

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network on which to apply the calculation
    sigma : dict key (string)
        The dictionary key containing the surface tension values to be used. If
        a pore property is given, it is interpolated to a throat list.
    theta : dict key (string)
        The dictionary key containing the contact angle values to be used. If
        a pore property is given, it is interpolated to a throat list.
    throat_diameter : dict key (string)
        The dictionary key containing the throat diameter values to be used.
    r_toroid : float or array_like
        The radius of the toroid surrounding the pore

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

    .. [1] G. Mason, N. R. Morrow, Effect of contact angle on capillary displacement
           curvatures in pore throats formed by spheres. J. Colloid Interface
           Sci. 168, 130 (1994).
    .. [2] J. Gostick, Random pore network modeling of fibrous PEMFC gas diffusion
           media using Voronoi and Delaunay tessellations. J. Electrochem.
           Soc. 160, F731 (2013).

    TODO: Triple check the accuracy of this equation
    """

    if surface_tension.split('.')[0] == 'pore':
        sigma = phase[surface_tension]
        sigma = phase.interpolate_data(data=sigma)
    else:
        sigma = phase[surface_tension]
    if contact_angle.split('.')[0] == 'pore':
        theta = phase[contact_angle]
        theta = phase.interpolate_data(data=theta)
    else:
        theta = phase[contact_angle]
    r = network[throat_diameter]/2
    R = r_toroid
    alpha = theta - 180 + _sp.arcsin(_sp.sin(_sp.radians(theta)/(1+r/R)))
    value = (-2*sigma/r) * \
        (_sp.cos(_sp.radians(theta - alpha)) /
            (1 + R/r*(1-_sp.cos(_sp.radians(alpha)))))
    value = value[phase.throats(physics.name)]
    return value
