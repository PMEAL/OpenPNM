r"""

.. autofunction:: openpnm.models.physics.meniscus.sinusoidal
.. autofunction:: openpnm.models.physics.meniscus.purcell
.. autofunction:: openpnm.models.physics.meniscus.general_toroidal

"""

import numpy as np
import logging
import sympy as syp
from openpnm.models.physics.capillary_pressure import _get_key_props
logger = logging.getLogger(__name__)


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
    network = target.project.network
    Dp_av = np.mean(network[pore_diameter][network['throat.conns']], axis=1)
    scale_b = r_toroid*network[throat_diameter]/Dp_av
    target['throat.scale_a'] = r_toroid
    target['throat.scale_b'] = scale_b
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

    throat_scale_a : dict key (string)
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
    x, a, b, rt, sigma, theta = syp.symbols('x, a, b, rt, sigma, theta')
    if profile_equation == 'elliptical':
        y = syp.sqrt(1 - (x/a)**2)*b
    elif profile_equation == 'sinusoidal':
        y = (syp.cos(syp.pi*x/(2*a)))*b
    else:
        logger.error('Profile equation is not valid')
    # Throat radius profile
    r = rt + (b-y)
    # Derivative of profile
    rprime = r.diff(x)
    # Filling angle
    alpha = syp.atan(rprime)
    # Radius of curvature of meniscus
    rm = r/syp.cos(alpha+theta)
    # distance from center of curvature to meniscus contact point (Pythagoras)
    d = syp.sqrt(rm**2 - r**2)
    # angle between throat axis, meniscus center and meniscus contact point
    gamma = syp.atan(r/d)
    # Capillary Pressure
    p = -2*sigma*syp.cos(alpha+theta)/r
    # Callable functions
    rx = syp.lambdify((x, a, b, rt), r, 'numpy')
    fill_angle = syp.lambdify((x, a, b, rt), alpha, 'numpy')
    Pc = syp.lambdify((x, a, b, rt, sigma, theta), p, 'numpy')
    rad_curve = syp.lambdify((x, a, b, rt, sigma, theta), rm, 'numpy')
    c2x = syp.lambdify((x, a, b, rt, sigma, theta), d, 'numpy')
    cap_angle = syp.lambdify((x, a, b, rt, sigma, theta), gamma, 'numpy')
    # All relative positions along throat
#    pos = np.arange(-0.999, 0.999, 1/num_points)
    hp = int(num_points/2)
    log_pos = np.logspace(-4, -1, hp+1)[:-1]
    lin_pos = np.arange(0.1, 1.0, 1/hp)
    half_pos = np.concatenate((log_pos, lin_pos))
    pos = np.concatenate((-half_pos[::-1], half_pos))
    # Now find the positions of the menisci along each throat axis
    Y, X = np.meshgrid(throatRad, pos)
    X *= fa
    # throat Capillary Pressure
    t_Pc = Pc(X, fa, fb, Y, surface_tension, contact)
    # Values of minima and maxima
    Pc_min = np.min(t_Pc, axis=0)
    Pc_max = np.max(t_Pc, axis=0)
    # Arguments of minima and maxima
    a_min = np.argmin(t_Pc, axis=0)
    a_max = np.argmax(t_Pc, axis=0)
    if mode == 'max':
        return Pc_max
    elif mode == 'touch':
        all_rad = rad_curve(X, fa, fb, Y, surface_tension, contact)
        all_c2x = c2x(X, fa, fb, Y, surface_tension, contact)
        all_cen = X + np.sign(all_rad)*all_c2x
        dist = all_cen + np.abs(all_rad)
        # Only count lengths where meniscus bulges into pore
        dist[all_rad > 0] = 0.0
        touch_len = network[touch_length]
        mask = dist > touch_len
        arg_touch = np.argmax(mask, axis=0)
        # Make sure we only count ones that happen before max pressure
        # And above min pressure (which will be erroneous)
        arg_in_range = (arg_touch < a_max) * (arg_touch > a_min)
        arg_touch[~arg_in_range] = a_max[~arg_in_range]
        x_touch = pos[arg_touch]*fa
        # Return the pressure at which a touch happens
        Pc_touch = Pc(x_touch, fa, fb, throatRad, surface_tension, contact)
        return Pc_touch
    elif target_Pc is None:
        logger.exception(msg='Please supply a target capillary pressure' +
                         ' when mode is "men"')
    if np.abs(target_Pc) < 1.0:
        target_Pc = 1.0
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
    men_data['c2x'] = c2x(xpos, fa, fb, throatRad, surface_tension, contact)
    men_data['gamma'] = cap_angle(xpos, fa, fb, throatRad,
                                  surface_tension, contact)
    men_data['radius'] = rad_curve(xpos, fa, fb, throatRad,
                                   surface_tension, contact)
    # xpos is relative to the throat center
    men_data['center'] = (xpos + np.sign(men_data['radius'])*men_data['c2x'])
    men_data['men_max'] = men_data['center'] - men_data['radius']
    logger.info(mode+' calculated for Pc: '+str(target_Pc))
    return men_data
