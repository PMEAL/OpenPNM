r"""
===============================================================================
Submodule -- capillary_pressure
===============================================================================

"""

import scipy as _sp


def washburn(target, surface_tension='pore.surface_tension',
             contact_angle='pore.contact_angle',
             diameter='throat.diameter'):
    r"""
    Computes the capillary entry pressure assuming the throat in a cylindrical
    tube.

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    sigma : string
        The dictionary key containing the surface tension values to be used. If
        a pore property is given, it is interpolated to a throat list.

    theta : string
        The dictionary key containing the contact angle values to be used. If
        a pore property is given, it is interpolated to a throat list.

    diameter : string
        The dictionary key containing the throat diameter values to be used.

    Notes
    -----
    The Washburn equation is:

    .. math::
        P_c = -\frac{2\sigma(cos(\theta))}{r}

    This is the most basic approach to calculating entry pressure and is
    suitable for highly non-wetting invading phases in most materials.

    """
    network = target.project.network
    phase = target.project.find_phase(target)
    if surface_tension.split('.')[0] == 'pore':
        sigma = phase.interpolate_data(propname=surface_tension)
    else:
        sigma = phase[surface_tension]
    if contact_angle.split('.')[0] == 'pore':
        theta = phase.interpolate_data(propname=contact_angle)
    else:
        theta = phase[contact_angle]
    r = network[diameter]/2
    value = -2*sigma*_sp.cos(_sp.radians(theta))/r
    if diameter.split('.')[0] == 'throat':
        value = value[phase.throats(target.name)]
    else:
        value = value[phase.pores(target.name)]
    value[_sp.absolute(value) == _sp.inf] = 0
    return value


def purcell(target, r_toroid, surface_tension='pore.surface_tension',
            contact_angle='pore.contact_angle',
            diameter='throat.diameter'):
    r"""
    Computes the throat capillary entry pressure assuming the throat is a
    toroid.

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    r_toroid : float or array_like
        The radius of the toroid surrounding the pore

    sigma : dict key (string)
        The dictionary key containing the surface tension values to be used.
        If a pore property is given, it is interpolated to a throat list.

    theta : dict key (string)
        The dictionary key containing the contact angle values to be used.
        If a pore property is given, it is interpolated to a throat list.

    diameter : dict key (string)
        The dictionary key containing the throat diameter values to be used.

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

    .. [1] G. Mason, N. R. Morrow, Effect of contact angle on capillary
           displacement curvatures in pore throats formed by spheres. J.
           Colloid Interface Sci. 168, 130 (1994).
    .. [2] J. Gostick, Random pore network modeling of fibrous PEMFC gas
           diffusion media using Voronoi and Delaunay tessellations. J.
           Electrochem. Soc. 160, F731 (2013).

    """
    network = target.project.network
    phase = target.project.find_phase(target)
    if surface_tension.split('.')[0] == 'pore':
        sigma = phase[surface_tension]
        sigma = phase.interpolate_data(propname=surface_tension)
    else:
        sigma = phase[surface_tension]
    if contact_angle.split('.')[0] == 'pore':
        theta = phase[contact_angle]
        theta = phase.interpolate_data(propname=contact_angle)
    else:
        theta = phase[contact_angle]
    r = network[diameter]/2
    R = r_toroid
    alpha = theta - 180 + _sp.arcsin(_sp.sin(_sp.radians(theta)/(1+r/R)))
    value = (-2*sigma/r) * \
        (_sp.cos(_sp.radians(theta - alpha)) /
            (1 + R/r*(1 - _sp.cos(_sp.radians(alpha)))))
    if diameter.split('.')[0] == 'throat':
        value = value[phase.throats(target.name)]
    else:
        value = value[phase.pores(target.name)]
    return value
