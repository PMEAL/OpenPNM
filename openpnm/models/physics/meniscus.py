r"""
Pore-scale models related to the meniscus calculations in pore/throats.
"""
import logging
import numpy as np
from openpnm.models.physics.capillary_pressure import _get_key_props
logger = logging.getLogger(__name__)

__all__ = ["purcell", "sinusoidal", "general_toroidal"]


def purcell(target,
            mode='max',
            target_Pc=None,
            num_points=1e3,
            r_toroid=5e-6,
            throat_diameter='throat.diameter',
            touch_length='throat.touch_length',
            surface_tension='pore.surface_tension',
            contact_angle='pore.contact_angle'):
    r'''
    Wrapper for the general toroidal model to implement the Purcell equation
    for a torus with cylindrical profile and toroidal radius r_toroid.

    Notes
    -----
    This approach accounts for the converging-diverging nature of many throat
    types. Advancing the meniscus beyond the apex of the toroid requires an
    increase in capillary pressure beyond that for a cylindical tube of the
    same radius. The details of this equation are described by Mason and
    Morrow [1]_, and explored by Gostick [2]_ in the context of a pore network
    model.
    '''
    target['throat.scale_a'] = r_toroid
    target['throat.scale_b'] = r_toroid
    output = general_toroidal(target=target,
                              mode=mode,
                              profile_equation='elliptical',
                              target_Pc=target_Pc,
                              num_points=num_points,
                              throat_diameter=throat_diameter,
                              touch_length=touch_length,
                              surface_tension=surface_tension,
                              contact_angle=contact_angle)
    return output


def sinusoidal(target,
               mode='max',
               target_Pc=None,
               num_points=1e3,
               r_toroid=5e-6,
               throat_diameter='throat.diameter',
               pore_diameter='pore.diameter',
               touch_length='throat.touch_length',
               surface_tension='pore.surface_tension',
               contact_angle='pore.contact_angle'):
    r'''
    Wrapper for the general toroidal model to implement the a sinusoidal
    profile. The quarter-wavelength is equal to toroidal radius r_toroid.
    The amplitude is equal to the toroidal radius multiplied by the ratio of
    the throat radius and average connecting pore radius.

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
    '''
    target['throat.scale_a'] = r_toroid
    target['throat.scale_b'] = r_toroid
    output = general_toroidal(target=target,
                              mode=mode,
                              profile_equation='sinusoidal',
                              target_Pc=target_Pc,
                              num_points=num_points,
                              throat_diameter=throat_diameter,
                              touch_length=touch_length,
                              surface_tension=surface_tension,
                              contact_angle=contact_angle)
    return output


def general_toroidal(target,
                     profile_equation='elliptical',
                     mode='max',
                     target_Pc=None,
                     num_points=1e3,
                     throat_scale_a='throat.scale_a',
                     throat_scale_b='throat.scale_b',
                     throat_diameter='throat.diameter',
                     touch_length='throat.touch_length',
                     surface_tension='pore.surface_tension',
                     contact_angle='pore.contact_angle'):
    r'''
    The general model for meniscus properties inside a toroidal throat

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    profile_equation : string (Default is 'elliptical')
        options elliptical, sinusoidal.

    mode : string (Default is 'max')
        Determines what information to send back. Options are:
        'max' : the maximum capillary pressure along the throat axis
        'touch' : the maximum capillary pressure a meniscus can sustain before
                  touching a solid feature
        None : return the meniscus info for a target pressure

    target_Pc : float
        The target capillary pressure to return data for when mode is 'men'

    num_points : float (Default 100)
        The number of divisions to make along the profile length to assess the
        meniscus properties in order to find target pressures, touch lengths,
        minima and maxima.

    throat_scale_a : dict key (string)
        The dictionary key containing the scale factor for adjusting the
        profile along the throat axis (x).

    throat_scale_b : dict key (string)
        The dictionary key containing the scale factor for adjusting the
        profile perpendicular to the throat axis (y).

    throat_diameter : dict key (string)
        The dictionary key containing the throat diameter values to be used.

    touch_length : dict key (string)
        The dictionary key containing the maximum length that a meniscus can
        protrude into the connecting pore before touching a solid feature and
        therfore invading

    surface_tension : dict key (string)
        The dictionary key containing the surface tension values to be used. If
        a pore property is given, it is interpolated to a throat list.

    contact_angle : dict key (string)
        The dictionary key containing the contact angle values to be used. If
        a pore property is given, it is interpolated to a throat list.
    '''
    from sympy import symbols, lambdify
    from sympy import atan as sym_atan
    from sympy import cos as sym_cos
    from sympy import sin as sym_sin
    from sympy import sqrt as sym_sqrt
    from sympy import pi as sym_pi

    # Get data from dictionary keys
    network = target.project.network
    phase = target.project.find_phase(target)
    (element,
     surface_tension,
     contact) = _get_key_props(phase=phase,
                               diameter=throat_diameter,
                               surface_tension=surface_tension,
                               contact_angle=contact_angle)
    # Contact Angle in radians
    contact = np.deg2rad(contact)
    # Network properties
    throatRad = network[throat_diameter]/2
    # Scaling parameters for throat profile
    fa = target[throat_scale_a]
    fb = target[throat_scale_b]
    # Governing equations
    x, a, b, rt, sigma, theta = symbols('x, a, b, rt, sigma, theta')
    if profile_equation == 'elliptical':
        y = sym_sqrt(1 - (x/a)**2)*b
    elif profile_equation == 'sinusoidal':
        y = (sym_cos((sym_pi/2)*(x/a)))*b
    else:
        logger.error('Profile equation is not valid, default to elliptical')
        y = sym_sqrt(1 - (x/a)**2)*b
    # Throat radius profile
    r = rt + (b-y)
    # Derivative of profile
    rprime = r.diff(x)
    # Filling angle
    alpha = sym_atan(rprime)
    # Angle between y axis and contact point to meniscus center
    eta = sym_pi - alpha - theta
    gamma = sym_pi/2 - eta
    # Radius of curvature of meniscus
    rm = r/sym_cos(eta)
    # distance from center of curvature to meniscus contact point (Pythagoras)
    d = rm*sym_sin(eta)
    # angle between throat axis, meniscus center and meniscus contact point
    # Capillary Pressure
    p = 2*sigma/rm
    # Callable functions
    rx = lambdify((x, a, b, rt), r, 'numpy')
    fill_angle = lambdify((x, a, b, rt), alpha, 'numpy')
    rad_curve = lambdify((x, a, b, rt, theta), rm, 'numpy')
    c2x = lambdify((x, a, b, rt, theta), d, 'numpy')
    cap_angle = lambdify((x, a, b, rt, theta), gamma, 'numpy')
    Pc = lambdify((x, a, b, rt, theta, sigma), p, 'numpy')
    # All relative positions along throat
    hp = int(num_points/2)
    log_pos = np.logspace(-4, -1, hp+1)[:-1]
    lin_pos = np.arange(0.1, 1.0, 1/hp)
    half_pos = np.concatenate((log_pos, lin_pos))
    pos = np.concatenate((-half_pos[::-1], half_pos))
    # Now find the positions of the menisci along each throat axis
    Y, X = np.meshgrid(throatRad, pos)
    X *= fa
    # throat Capillary Pressure
    t_Pc = Pc(X, fa, fb, Y, contact, surface_tension)
    # Values of minima and maxima
    Pc_min = np.min(t_Pc, axis=0)
    Pc_max = np.max(t_Pc, axis=0)
    # Arguments of minima and maxima
    a_min = np.argmin(t_Pc, axis=0)
    a_max = np.argmax(t_Pc, axis=0)
    if mode == 'max':
        return Pc_max
    elif mode == 'touch':
        all_rad = rad_curve(X, fa, fb, Y, contact)
        all_c2x = c2x(X, fa, fb, Y, contact)
        all_cen = X - all_c2x
        dist = all_cen + all_rad
        # Only count lengths where meniscus bulges into pore
        dist[all_rad < 0] = 0.0
        touch_len = network[touch_length]
        mask = dist > touch_len
        arg_touch = np.argmax(mask, axis=0)
        # Make sure we only count ones that happen before max pressure
        # And above min pressure (which will be erroneous)
        arg_in_range = (arg_touch < a_max) * (arg_touch > a_min)
        arg_touch[~arg_in_range] = a_max[~arg_in_range]
        x_touch = pos[arg_touch]*fa
        # Return the pressure at which a touch happens
        Pc_touch = Pc(x_touch, fa, fb, throatRad, contact, surface_tension)
        return Pc_touch
    elif target_Pc is None:
        logger.error(msg='Please supply a target capillary pressure'
                     + ' when mode is "men", default to 1.0e-6')
        target_Pc = 1.0e-6
    if np.abs(target_Pc) < 1.0e-6:
        logger.error(msg='Please supply a target capillary pressure'
                     + ' with absolute value greater than 1.0e-6,'
                     + ' default to 1.0e-6')
        target_Pc = 1.0e-6
    # Find the position in-between the minima and maxima corresponding to
    # the target pressure
    inds = np.indices(np.shape(t_Pc))
    # Change values outside the range between minima and maxima to be those
    # Values
    mask = inds[0] < np.ones(len(pos))[:, np.newaxis]*a_min
    t_Pc[mask] = (np.ones(len(pos))[:, np.newaxis]*Pc_min)[mask]
    mask = inds[0] > np.ones(len(pos))[:, np.newaxis]*a_max
    t_Pc[mask] = (np.ones(len(pos))[:, np.newaxis]*Pc_max)[mask]
    # Find the argument at or above the target Pressure
    mask = t_Pc >= target_Pc
    arg_x = np.argmax(mask, axis=0)
    # If outside range change to minima or maxima accordingly
    arg_x[target_Pc < Pc_min] = a_min[target_Pc < Pc_min]
    arg_x[target_Pc > Pc_max] = a_max[target_Pc > Pc_max]
    xpos = pos[arg_x]*fa
    xmin = pos[a_min]*fa
    xmax = pos[a_max]*fa
    # Output
    men_data = {}
    men_data['pos'] = xpos
    men_data['rx'] = rx(xpos, fa, fb, throatRad)
    men_data['alpha'] = fill_angle(xpos, fa, fb, throatRad)
    men_data['alpha_min'] = fill_angle(xmin, fa, fb, throatRad)
    men_data['alpha_max'] = fill_angle(xmax, fa, fb, throatRad)
    men_data['c2x'] = c2x(xpos, fa, fb, throatRad, contact)
    men_data['gamma'] = cap_angle(xpos, fa, fb, throatRad, contact)
    men_data['radius'] = rad_curve(xpos, fa, fb, throatRad, contact)
    # xpos is relative to the throat center
    men_data['center'] = (xpos - men_data['c2x'])
    men_data['men_max'] = men_data['center'] - men_data['radius']

    logger.info(mode+' calculated for Pc: '+str(target_Pc))
    return men_data
