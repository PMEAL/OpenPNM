from scipy import ndimage
import scipy as _sp
import numpy as np
import logging
import sympy as syp
import pandas as pd
from openpnm.models.physics.capillary_pressure import _get_key_props
logger = logging.getLogger(__name__)


def _prop_parser(obj, prop, entity):
    r'''
    Helper function to get data in pore or throat format depending on what
    you want
    '''
    if (prop.split('.')[0] == 'pore' and
       entity.split('.')[0] == 'throat'):
        value = obj[prop]
        value = obj.interpolate_data(propname=prop)
    else:
        value = obj[prop]
    return value


def toroidal(target,
             mode='max',
             r_toroid=5e-6,
             target_Pc=None,
             surface_tension='pore.surface_tension',
             contact_angle='pore.contact_angle',
             diameter='throat.diameter'):
    r"""
    Calculate the filling angle (alpha) for a given capillary pressure

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.
    mode : string (Default is 'max')
        Determines what information to send back. Options are:
        'max' : the maximum capillary pressure along the throat axis
        'men' : return the meniscus info for a target pressure
    r_toroid : float or array_like
        The radius of the toroid surrounding the pore
    target_Pc : float
        The target capillary pressure
    surface_tension : dict key (string)
        The dictionary key containing the surface tension values to be used. If
        a pore property is given, it is interpolated to a throat list.
    contact_angle : dict key (string)
        The dictionary key containing the contact angle values to be used. If
        a pore property is given, it is interpolated to a throat list.
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

    !!! Takes mean contact angle and surface tension !!!

    """
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
        check_max = fill_angle > a_max
        if np.any(check_max):
            if np.size(a_max) > 1:
                fill_angle[check_max] = a_max[check_max]
            else:
                fill_angle[check_max] = a_max
        r_men = R*(1+(ratio)-_sp.cos(fill_angle))/_sp.cos(theta-fill_angle)
        Pc = 2*sigma/r_men
        return Pc

    Pc_max = purcell_pressure(ratios, a_max, theta, sigma, R)
    if mode == 'max':
        return Pc_max
    elif target_Pc is None:
        logger.exception(msg='Please supply a target capillary pressure')
    else:
        pass
    fill_angle = _sp.deg2rad(np.linspace(-30, 150, 1001))

    alpha = np.zeros_like(ratios)
    for T, ratio in enumerate(ratios):
        mask = np.zeros_like(fill_angle, dtype=bool)
        nudge = 100
        all_Pc = purcell_pressure(ratio, fill_angle, theta[T], sigma[T], R)
        if target_Pc > all_Pc.max():
            # Target Pc out of range
            lowest = fill_angle[np.argwhere(all_Pc == all_Pc.max())[0][0]]
        else:
            while np.sum(mask) == 0:
                plus_mask = all_Pc < target_Pc + nudge
                minus_mask = all_Pc > target_Pc - nudge
                mask = np.logical_and(plus_mask, minus_mask)
                if np.sum(mask) == 0:
                    nudge += 10

            regions = ndimage.find_objects(ndimage.label(mask)[0])
            rx = [np.mean(fill_angle[regions[r]]) for r in range(len(regions))]
            root_x = np.asarray(rx)
            lowest = np.min(root_x)
        alpha[T] = lowest

    logger.info('Filling angles calculated for Pc: '+str(target_Pc))
    men_data = {}
    f = theta-alpha
    # Handle potential divide by zero
    f[np.abs(f) == np.pi/2] = f[np.abs(f) == np.pi/2]*(1-1e-12)
    # Meniscus radius
    r_men = R*(1 + rt/R - np.cos(alpha)) / np.cos((f))
    # Vertical adjustment for centre of circle
    y_off = r_toroid*np.sin(alpha)
    # Angle between contact point - centre - vertical
    zeta = (theta-alpha-np.pi/2)
    # Distance that center of meniscus is below the plane of the throat
    center = y_off - r_men*np.cos(zeta)
    men_data['alpha'] = _sp.rad2deg(alpha)
    men_data['alpha_max'] = a_max
    men_data['radius'] = r_men
    men_data['center'] = center
    return men_data


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

    def in_range(target_Pc):
        r'''
        Check whether the target pressure is in range for each throat
        '''
        return np.logical_and((target_Pc >= min_Pc), (target_Pc <= max_Pc))

    def get_root(target_Pc):
        r'''
        Get the root between the minima and maxima
        '''
        # interpolated initial guess
        x0 = min_point+(x_range*(target_Pc-min_Pc)/(Pc_range))
        x0[~in_range(target_Pc)] = np.nan
        # find root with function adjusted for target
        root = Newton_Raphson(x0,
                              poreRad,
                              throatRad,
                              throatLength,
                              sigma,
                              np.deg2rad(theta),
                              target_Pc)
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
    men_data['radius'] = rad_curve(pos, poreRad, throatRad, throatLength,
                                   sigma, theta, offset)
    men_data['center'] = pos - np.sign(target_Pc)*men_data['alpha']
    logger.info(mode+' calculated for Pc: '+str(target_Pc))
    return men_data
