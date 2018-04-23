r"""
===============================================================================
Submodule -- capillary_pressure
===============================================================================

"""

import scipy as _sp
import sympy as syp
import numpy as np
import pandas as pd
import logging
from transforms3d import _gohlketransforms as tr
logger = logging.getLogger(__name__)


def _get_key_props(phase=None, diameter='throat.diameter',
                   surface_tension='pore.surface_tension',
                   contact_angle='pore.contact_angle'):
    r"""
    Many of the methods are generic to pores and throats. Some information may
    be stored on either the pore or throat and needs to be interpolated.
    This is a helper method to return the properties in the correct format.
    To do:
        Check for method to convert throat to pore data
    """
    element = diameter.split('.')[0]
    if element == 'pore':
        if 'throat' in surface_tension:
            sigma = phase.interpolate_data(propname=surface_tension)
        else:
            sigma = phase[surface_tension]
        if 'throat' in contact_angle:
            theta = phase.interpolate_data(propname=contact_angle)
        else:
            theta = phase[contact_angle]
    if element == 'throat':
        if 'pore' in surface_tension:
            sigma = phase.interpolate_data(propname=surface_tension)
        else:
            sigma = phase[surface_tension]
        if 'pore' in contact_angle:
            theta = phase.interpolate_data(propname=contact_angle)
        else:
            theta = phase[contact_angle]

    return element, sigma, theta


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
    element, sigma, theta = _get_key_props(phase=phase,
                                           diameter=diameter,
                                           surface_tension=surface_tension,
                                           contact_angle=contact_angle)
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
    element, sigma, theta = _get_key_props(phase=phase,
                                           diameter=diameter,
                                           surface_tension=surface_tension,
                                           contact_angle=contact_angle)
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


def sinusoidal(target,
               mode='max',
               target_Pc=None,
               surface_tension='pore.surface_tension',
               contact_angle='pore.contact_angle',
               pore_diameter='pore.diameter',
               throat_diameter='throat.diameter',
               throat_length='throat.length',
               **kwargs):
    r"""
    The profile of a throat is approximated with a sinusoidal function
    that depends on the average of the connecting pore diameters and throat
    diamater. It represents a converging-diverging geometry that has a minima
    at the mid-point of the throat and produces similar behaviour to the
    Purcell model but allows for a more slowly varying profile at the
    ends of the throat.

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.
    mode : string (Default is 'max')
        Determines what information to send back. Options are:
        'max' : the maximum capillary pressure along the throat axis, does not
        'men' : return the meniscus info for a target pressure
    target_Pc : float (Default is None)
        The target capillary pressure for use with mode 'men'
    surface_tension : dict key (string)
        The dictionary key containing the surface tension values to be used. If
        a pore property is given, it is interpolated to a throat list.
    contact_angle : dict key (string)
        The dictionary key containing the contact angle values to be used. If
        a pore property is given, it is interpolated to a throat list.
    pore_diameter : dict key (string)
        The dictionary key containing the pore diameter values to be used.
    throat_diameter : dict key (string)
        The dictionary key containing the throat diameter values to be used.
    throat_length : dict key (string)
        The dictionary key containing the throat length values to be used.

    Notes
    -----
    The capillary pressure equation for a sinusoidal throat is extended from
    the Washburn equation as [1]_:

    .. math::
        P_c = -\frac{2\sigma(cos(\alpha + \theta))}{r(x)}

    where alpha is:
    .. math::
        \alpha = arctan(\frac{dr}{dx})

    References
    ----------

    .. [1] A. Forner-Cuenca et. al, Advanced Water Management in PEFCs:
        Diffusion Layers with Patterned Wettability.
        J. ECS. 163, 9, F1038-F1048 (2016).
    """
    network = target.project.network
    phase = target.project.find_phase(target)
    element, sigma, theta = _get_key_props(phase=phase,
                                           diameter=throat_diameter,
                                           surface_tension=surface_tension,
                                           contact_angle=contact_angle)
    # Symbols
    # sigma
    s = syp.Symbol('s')
    # theta
    t = syp.Symbol('t')
    # position of mensicus along throat axis, zero at center
    x = syp.Symbol('x')
    # pore radius
    rp = syp.Symbol('rp')
    # throat radius
    rt = syp.Symbol('rt')
    # throat lenggth
    l = syp.Symbol('l')
    # Pressure offset for finding minima
    off = syp.Symbol('off')
    # Equations
    # Radius profile along throat length
    y = (rp-rt)*(1-syp.cos(2*syp.pi*x/l))/2 + rt
    # dr/dx used for filling angle
    yprime = y.diff(x)
    # Filling angle
    alpha = syp.atan(yprime)
    # Meniscus Radius of curvature
    R = y/syp.cos(alpha+t)
    # distance from center of curvature to meniscus contact point (Pythagoras)
    a = syp.sqrt(R*R - y*y)
    # angle between throat axis, meniscus center and meniscus contact point
    gamma = syp.asin(y/R)
    # Capillary Pressure function with target capillary pressure adjustment for
    # root finding
    f = -2*s*syp.cos(alpha+t)/y - off
    # df/dx used for Newton-Raphson method for root finding
    fprime = f.diff(x)
    # Callable expressions
    rx = syp.lambdify((x, rp, rt, l), y, 'numpy')
    Pc = syp.lambdify((x, rp, rt, l, s, t, off), f, 'numpy')
    Pc_prime = syp.lambdify((x, rp, rt, l, s, t, off), fprime, 'numpy')
    rad_curve = syp.lambdify((x, rp, rt, l, s, t, off), R, 'numpy')
    c2x = syp.lambdify((x, rp, rt, l, s, t, off), a, 'numpy')
    fill_angle = syp.lambdify((x, rp, rt, l), alpha, 'numpy')
    cap_angle = syp.lambdify((x, rp, rt, l, s, t, off), gamma, 'numpy')
    # Network properties
    throatLength = network[throat_length]
    poreRad = np.mean(network[pore_diameter][network['throat.conns']], axis=1)
    poreRad /= 2
    throatRad = network[throat_diameter]/2
    Nt = network.Nt
    # Model ouputs
    offset = np.zeros(Nt)
    min_Pc = np.zeros(Nt)
    max_Pc = np.zeros(Nt)
    min_arg = np.zeros(Nt, dtype=int)
    max_arg = np.zeros(Nt, dtype=int)
    min_point = np.zeros(Nt)
    max_point = np.zeros(Nt)

    # Preprocessing - go along the throat length and work out min and max Pc
    # and the position where this occurs
    for i in range(Nt):
        points = np.arange(-throatLength[i]/2,
                           throatLength[i]/2,
                           throatLength[i]/100)
        all_Pc = Pc(points,
                    poreRad[i],
                    throatRad[i],
                    throatLength[i],
                    sigma[i],
                    np.deg2rad(theta[i]),
                    offset[i])
        min_Pc[i] = np.min(all_Pc)
        max_Pc[i] = np.max(all_Pc)
        min_arg[i] = np.argmin(all_Pc)
        max_arg[i] = np.argmax(all_Pc)
        min_point[i] = points[min_arg[i]]
        max_point[i] = points[max_arg[i]]
    if mode == 'max':
        return max_Pc
    elif target_Pc is None:
        logger.exception(msg='Please supply a target capillary pressure')
    else:
        pass
    # If we got this far we're looking for the meniscus information for a
    # target Capillary pressure
    Pc_range = max_Pc-min_Pc
    x_range = max_point - min_point

    # Private helper functions
    def Newton_Raphson(x0, rp, rt, l, s, t, off):
        tol = np.ones(len(x0), dtype=float)
        n = 0
        n_max = 50
        while np.any(tol[~np.isnan(tol)] > 1e-6) and n < n_max:
            func = Pc(x0, rp, rt, l, s, t, off)
            func_prime = Pc_prime(x0, rp, rt, l, s, t, off)
            xn = x0 - func/func_prime
            tol = np.abs((xn - x0)/x0)
            x0 = xn
            n += 1
        return x0

    def in_range(target):
        r'''
        Check whether the target pressure is in range for each throat
        '''
        return np.logical_and((target >= min_Pc), (target <= max_Pc))

    def get_root(target):
        r'''
        Get the root between the minima and maxima
        '''
        # interpolated initial guess
        x0 = min_point+(x_range*(target-min_Pc)/(Pc_range))
        x0[~in_range(target)] = np.nan
        # find root with function adjusted for target
        root = Newton_Raphson(x0,
                              poreRad,
                              throatRad,
                              throatLength,
                              sigma,
                              np.deg2rad(theta),
                              target)
        return root

    # Now find the positions of the menisci along each throat axis
    men_data = {}
    pos = get_root(target_Pc)
    men_data['pos'] = pos
    men_data['rx'] = rx(pos, poreRad, throatRad, throatLength)
    men_data['alpha'] = fill_angle(pos, poreRad, throatRad, throatLength)
    men_data['beta'] = c2x(pos, poreRad, throatRad, throatLength,
                           sigma, theta, offset)
    men_data['gamma'] = cap_angle(pos, poreRad, throatRad, throatLength,
                                  sigma, theta, offset)
    men_data['rad'] = rad_curve(pos, poreRad, throatRad, throatLength,
                                sigma, theta, offset)
    men_data['cen'] = pos - np.sign(target_Pc)*men_data['alpha']
    df = pd.DataFrame(men_data)
    rec_arr = df.to_records(index=False)
    logger.info(mode+' calculated for Pc: '+str(target_Pc))
    return rec_arr


def ransohoff_snap_off(target,
                       shape_factor=2.0,
                       require_pair=True,
                       contact_angle='pore.contact_angle',
                       surface_tension='pore.surface_tension',
                       diameter='throat.diameter',
                       wavelength=5e-6,
                       vertices='throat.offset_vertices',
                       **kwargs):
    r"""
    Computes the capillary snap-off pressure assuming the throat is cylindrical
    with converging-diverging change in diamater - like the Purcell model.
    The wavelength of the change in diamater is the fiber radius.
    Ref: Ransohoff, T.C., Gauglitz, P.A. and Radke, C.J., 1987. Snap‐off of gas
    bubbles in smoothly constricted noncircular capillaries. AIChE Journal,
    33(5), pp.753-765.

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.
    shape_factor :
        constant dependent on the shape of throat cross-section 1.75 - 2.0, see
        Ref
    sigma : dict key (string)
        The dictionary key containing the surface tension values to be used. If
        a pore property is given, it is interpolated to a throat list.
    throat_diameter : dict key (string)
        The dictionary key containing the throat diameter values to be used.
    wavelength :
        The dictionary key containing the radius of the transverse interfacial
        radius of curvature at the neck (fiber radius in fibrous media)

    Notes
    -----
    This equation should be used to calculate the snap off capillary pressure
    in fribrous media

    """
    network = target.project.network
    phase = target.project.find_phase(target)
    element, sigma, theta = _get_key_props(phase=phase,
                                           diameter=diameter,
                                           surface_tension=surface_tension,
                                           contact_angle=contact_angle)
    try:
        geometry = network.geometries(network.geometries()[0])[0]
        all_verts = geometry[vertices]
        # Work out whether throat geometry can support at least one pair of
        # adjacent arc menisci that can grow and merge to form snap-off
        # Only works if throat vertices are in convex hull order
        angles_ok = np.zeros(network.Nt, dtype=bool)
        for T in range(network.Nt):
            verts = all_verts[T]
            x = verts[:, 0]
            y = verts[:, 1]
            z = verts[:, 2]
            # PLus
            p = 1
            # Minus
            m = -1
            verts_p = np.vstack((np.roll(x, p),
                                 np.roll(y, p),
                                 np.roll(z, p))).T
            verts_m = np.vstack((np.roll(x, m),
                                 np.roll(y, m),
                                 np.roll(z, m))).T
            v1 = verts_p - verts
            v2 = verts_m - verts
            corner_angles = np.rad2deg(tr.angle_between_vectors(v1,
                                                                v2,
                                                                axis=1))
            # Logical test for existence of arc menisci
            am = theta[T] <= 90 - corner_angles/2
            if require_pair:
                # Logical test for two adjacent arc menisci
                pair_p = np.logical_and(am, np.roll(am, + p))
                pair_m = np.logical_and(am, np.roll(am, + m))
                am_pair = np.any(np.logical_or(pair_p, pair_m))
                angles_ok[T] = am_pair
            else:
                # Logical test for any arc menisci
                angles_ok[T] = np.any(am)
    except:
        logger.warning("Model is designed to work with property: " +
                       vertices)
        angles_ok = np.ones(network.Nt, dtype=bool)

    # Condition for arc menisci to form in corners
    rad_Ts = network[diameter]/2
    # Ransohoff and Radke eq. 4
    C = 1/rad_Ts - 1/wavelength
    value = sigma*C
    # Only throats that can support arc menisci can snap-off
    value[~angles_ok] = np.nan
    logger.info("Snap off pressures calculated for " +
                str(np.around(100*np.sum(angles_ok)/np.size(angles_ok), 0)) +
                "% of throats")
    return value[phase.throats(target.name)]


def purcell_bi(target, r_toroid,
               surface_tension='pore.surface_tension',
               contact_angle='pore.contact_angle',
               diameter='throat.diameter',
               h_max='pore.diameter',
               max_dist=True,
               **kwargs):
    r"""
    Computes the throat capillary entry pressure assuming the throat is a
    toroid. Pressure is a function of both pore and throat diameters so it
    becomes a bi-dirctional calculation.

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

    """
    network = target.project.network
    phase = target.project.find_phase(target)
    element, sigma, theta = _get_key_props(phase=phase,
                                           diameter=diameter,
                                           surface_tension=surface_tension,
                                           contact_angle=contact_angle)
    # Mason and Morrow have the definitions switched
    theta = 180 - theta
    th = _sp.deg2rad(theta)
    rt = network[diameter]/2
    R = r_toroid
    a_max = th - np.arcsin((np.sin(th))/(1+rt/R))
    if max_dist and element == 'throat':
        # Perform analysis for entry into both pores
        r_max = np.zeros([network.Nt, 2])
        for j in range(2):
            Pj = network['throat.conns'][:, j]
            dj = network[h_max][Pj]
            max_reached = np.zeros(network.Nt, dtype=bool)
            alpha_reached = np.zeros(network.Nt)
            # With increasing filling angle assess whether interface has passed
            # the critical distance and record the critical angle
            a_space = _sp.linspace(1e-3, _sp.pi, 181)
            nudge = 0.001
            for a_test in a_space:
                nudgers = network.throats()[np.around(th-a_test, 0) == 90]
                if len(nudgers) > 0:
                    th[nudgers] += nudge
                r = R*(1+(rt/R)-_sp.cos(a_test))/_sp.cos(th-a_test)
                # Vertical adjustment for centre of circle
                y_off = R*np.sin(a_test)
                # Angle between contact point - centre - vertical
                zeta = (th-a_test-(np.pi/2))
                c = y_off - r*np.cos(zeta)
                y_max = c+r
                ts = network.throats()[(y_max > (dj)) * (~max_reached)]
                if len(ts) > 0:
                    max_reached[ts] = True
                    alpha_reached[ts] = a_test
                if len(nudgers) > 0:
                    th[nudgers] -= nudge
            # Any interfaces that never reach a wall are ok"
            alpha_reached[~max_reached] = a_max[~max_reached]
            temp = max_reached[np.abs(a_max) > np.abs(alpha_reached)]
            temp_sum = np.around(100*np.sum(temp)/len(a_max), 2)
            logger.info("Percentage max before BT " + str(temp_sum))
            # Any interfaces that can expand to maximum curvature before
            # hitting a pore wall are ok
            mask = np.abs(a_max) < np.abs(alpha_reached)
            alpha_reached[mask] = a_max[mask]
            r_max[:, j] = (R*(1 + (rt/R) - _sp.cos(alpha_reached)) /
                           _sp.cos(th - alpha_reached))

        value = (2*np.vstack((sigma, sigma)).T) / r_max
    else:
        r_max = R*(1+(rt/R)-_sp.cos(a_max))/_sp.cos(th-a_max)
        value = 2*sigma/r_max

    if element == 'throat':
        value = value[phase.throats(target.name)]
    else:
        value = value[phase.pores(target.name)]
    return value


def purcell_filling_angle(target, r_toroid,
                          surface_tension='pore.surface_tension',
                          contact_angle='pore.contact_angle',
                          diameter='throat.diameter',
                          Pc=1e3,
                          **kwargs):
    r"""
    Calculate the filling angle (alpha) for a given capillary pressure

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

    !!! Takes mean contact angle and surface tension !!!

    """
    from scipy import ndimage
    network = target.project.network
    phase = target.project.find_phase(target)
    element, sigma, theta = _get_key_props(phase=phase,
                                           diameter=diameter,
                                           surface_tension=surface_tension,
                                           contact_angle=contact_angle)
    # Mason and Morrow have the definitions switched
    theta = 180 - theta
    theta = _sp.deg2rad(theta)
    rt = network[diameter]/2
    R = r_toroid
    ratios = rt/R
    a_max = theta - np.arcsin((np.sin(theta))/(1 + ratios))

    def purcell_pressure(ratio, fill_angle, theta, sigma, R):
        # Helper function
        a_max = theta - np.arcsin((np.sin(theta))/(1+ratio))
        fill_angle[fill_angle > a_max] = a_max
        r_men = R*(1+(ratio)-_sp.cos(fill_angle))/_sp.cos(theta-fill_angle)
        Pc = 2*sigma/r_men
        return Pc

    fill_angle = _sp.deg2rad(np.linspace(-30, 150, 1001))

    alpha = np.zeros_like(ratios)
    for T, ratio in enumerate(ratios):
        mask = np.zeros_like(fill_angle, dtype=bool)
        nudge = 100
        all_Pc = purcell_pressure(ratio, fill_angle, theta[T], sigma[T], R)
        if Pc > all_Pc.max():
            # Target Pc out of range
            lowest = fill_angle[np.argwhere(all_Pc == all_Pc.max())[0][0]]
        else:
            while np.sum(mask) == 0:
                plus_mask = all_Pc < Pc + nudge
                minus_mask = all_Pc > Pc - nudge
                mask = np.logical_and(plus_mask, minus_mask)
                if np.sum(mask) == 0:
                    nudge += 10

            regions = ndimage.find_objects(ndimage.label(mask)[0])
            rx = [np.mean(fill_angle[regions[r]]) for r in range(len(regions))]
            root_x = np.asarray(rx)
            lowest = np.min(root_x)
        alpha[T] = lowest

    logger.info('Filling angles calculated for Pc: '+str(Pc))
    target['throat.alpha_max'] = a_max
    return _sp.rad2deg(alpha)


def _prop_parser(obj, prop, entity):
    r'''
    Helper function to get data in pore or throat format depending on what
    you want
    '''
    if (prop.split('.')[0] == 'pore' and
       entity.split('.')[0] == 'throat'):
        value = obj.interpolate_data(prop)
    else:
        value = obj[prop]
    return value


def purcell_meniscus_radius(target, r_toroid,
                            contact_angle='pore.contact_angle',
                            filling_angle='throat.alpha',
                            diameter='throat.diameter',
                            **kwargs):
    r"""
    Function to return the radius of curvature for the sphere whose spherical
    cap forms the meniscus inside a throat as per the Purcell model.
    Assumes angles are stored in degrees
    """
    network = target.project.network
    phase = target.project.find_phase(target)
    theta = _prop_parser(phase, contact_angle, diameter)
    # Mason and Morrow have the definitions switched
    theta = 180 - theta
    alpha = _prop_parser(phase, filling_angle, diameter)
    theta = np.deg2rad(theta)
    alpha = np.deg2rad(alpha)
    R = r_toroid
    rt = network[diameter]/2
    f = theta-alpha
    # Handle potential divide by zero
    f[np.abs(f) == np.pi/2] = f[np.abs(f) == np.pi/2]*(1-1e-12)
    r_men = R*(1 + rt/R - np.cos(alpha)) / np.cos((f))
    return r_men


def purcell_meniscus_center(target, r_toroid,
                            contact_angle='pore.contact_angle',
                            filling_angle='throat.alpha',
                            men_rad='throat.meniscus_radius',
                            normal='throat.normal',
                            center='throat.centroid',
                            **kwargs):
    r"""
    Function to return the center offset of the sphere whose spherical
    cap forms the meniscus inside a throat as per the Purcell model.
    N.B to be multiplied by the throat normal vector
    """
    phase = target.project.find_phase(target)
    theta = _prop_parser(phase, contact_angle, men_rad)
    # Mason and Morrow have the definitions switched
    theta = 180 - theta
    alpha = _prop_parser(phase, filling_angle, men_rad)
    theta = np.deg2rad(theta)
    alpha = np.deg2rad(alpha)
    # Radius of meniscus between fibres
    r_men = target[men_rad]
    # Vertical adjustment for centre of circle
    y_off = r_toroid*np.sin(alpha)
    # Angle between contact point - centre - vertical
    zeta = (theta-alpha-np.pi/2)
    # Store this for coop filling analysis
    target['throat.zeta'] = zeta
    # Distance that center of meniscus is below the plane of the throat
    value = y_off - r_men*np.cos(zeta)

    return value
