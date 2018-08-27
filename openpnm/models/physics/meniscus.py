import scipy as _sp
import numpy as np
import logging
import sympy as syp
from openpnm.models.physics.capillary_pressure import _get_key_props
logger = logging.getLogger(__name__)


def toroidal(target,
             mode='max',
             r_toroid=5e-6,
             target_Pc=None,
             surface_tension='pore.surface_tension',
             contact_angle='pore.contact_angle',
             diameter='throat.diameter',
             touch_length='throat.touch_length'):
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
    touch_length : dict key (string)
        The dictionary key containing the maximum length that a meniscus can
        protrude into the connecting pore before touching a solid feature and
        therfore invading

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
    theta = _sp.deg2rad(theta)
    rt = network[diameter]/2
    Rf = r_toroid

    a, r, R, s, t, p = syp.symbols('a, r, R, s, t, p')
    rhs = -2*s/(R*(1+(r/R)-syp.cos(a))/syp.cos(t-a))
#    lhs = p
#    roots = syp.solve(lhs-rhs, a)
#    _, r1, r2, _ = roots
    eit = syp.exp(syp.I*t)
    # There are 2 non-trivial solutions when rearranging the Purcell Pc eqn
    # To solve for alpha (a), r1 and r2 - r2 is outside a_min and a_max,
    # r1 is inside the range. It takes sympy to solve the equations so the
    # root inside the range is provided below but can be verified by running
    # the commented out lines above
    r1 = -syp.I*syp.log((R*p*eit + p*r*eit -
                         syp.sqrt((2*R*p**2*r*eit +
                                   2*R*p*s*syp.exp(2*syp.I*t) +
                                   2*R*p*s + p**2*r**2*eit -
                                   4*s**2*eit)*eit))/(R*p*eit - 2*s))
    r2 = -syp.I*syp.log((R*p*eit + p*r*eit +
                         syp.sqrt((2*R*p**2*r*eit +
                                   2*R*p*s*syp.exp(2*syp.I*t) +
                                   2*R*p*s + p**2*r**2*eit -
                                   4*s**2*eit)*eit))/(R*p*eit - 2*s))
    a_min = t - syp.asin((syp.sin(t))/(1+r/R))
    a_max = t - syp.pi + syp.asin((syp.sin(t))/(1+r/R))
    # alpha at given Pc
    fa_Pc = syp.lambdify((p, r, R, s, t), r1, 'numpy')
    # Pc at given alpha
    fPc = syp.lambdify((a, r, R, s, t), rhs, 'numpy')
    # alphas where max and min Pc occurs
    fa_max = syp.lambdify((r, R, t), a_max, 'numpy')
    fa_min = syp.lambdify((r, R, t), a_min, 'numpy')
    # Values at min and max
    a_maxs = fa_max(rt, Rf, theta)
    pc_max = fPc(a_maxs, rt, Rf, sigma, theta)
    a_mins = fa_min(rt, Rf, theta)
    pc_min = fPc(a_mins, rt, Rf, sigma, theta)

    if mode == 'max':
        return pc_max
    elif mode == 'touch':
        pos = np.linspace(-np.pi/2, np.pi/2, 91)[:: -1]
        Y, X = np.meshgrid(rt, pos)
        f = theta-X
        # Handle potential divide by zero
        f[np.abs(f) == np.pi/2] = f[np.abs(f) == np.pi/2]*(1-1e-12)
        # Meniscus radius
        r_men = Rf*(1 + Y/Rf - np.cos(X)) / np.cos((f))
        # Vertical adjustment for centre of circle
        y_off = Rf*np.sin(X)
        # Angle between contact point - centre - vertical
        zeta = (f-np.pi/2)
        # Distance that center of meniscus is below the plane of the throat
        center = y_off - r_men*np.cos(zeta)
        dist = center + r_men
        # Only count lengths where meniscus bulges into pore
        dist[r_men > 0] = 0.0
        touch_len = network[touch_length]
        mask = dist < -touch_len
        arg_touch = np.argmax(mask, axis=0)
        # Make sure we only count ones that happen before max pressure
        # And above min pressure (which will be erroneous)
        x_touch = pos[arg_touch]
        arg_in_range = (x_touch > a_maxs) * (x_touch < a_mins)
        x_touch[~arg_in_range] = x_touch[~arg_in_range]
        # Return the pressure at which a touch happens
        Pc_touch = fPc(x_touch, rt, Rf, sigma, theta)
        return Pc_touch
    elif target_Pc is None:
        logger.exception(msg='Please supply a target capillary pressure' +
                         ' when mode is "men"')
    if np.abs(target_Pc) < 1.0:
        target_Pc = 1.0
    # Masks to determine which throats to actually calculate alpha for
    # Outside the valid range of pressures min or max values are used
    over_range = target_Pc > pc_max
    undr_range = target_Pc < pc_min
    in_range = ~over_range * ~undr_range
    alpha = np.zeros(len(rt))
    if np.any(in_range):
        alpha[in_range] = np.real(fa_Pc(target_Pc, rt[in_range], Rf,
                                        sigma[in_range], theta[in_range]))
    if np.any(over_range):
        alpha[over_range] = a_maxs[over_range]
    if np.any(undr_range):
        alpha[undr_range] = a_mins[undr_range]

    logger.info('Filling angles calculated for Pc: '+str(target_Pc))
    men_data = {}
    f = theta-alpha
    # Handle potential divide by zero
    f[np.abs(f) == np.pi/2] = f[np.abs(f) == np.pi/2]*(1-1e-12)
    # Meniscus radius
    r_men = Rf*(1 + rt/Rf - np.cos(alpha)) / np.cos((f))
    # Vertical adjustment for centre of circle
    y_off = Rf*np.sin(alpha)
    # Angle between contact point - centre - vertical
    zeta = (theta-alpha-np.pi/2)
    # Distance that center of meniscus is below the plane of the throat
    center = y_off - r_men*np.cos(zeta)
    men_data['alpha'] = alpha
    men_data['alpha_max'] = a_maxs
    men_data['alpha_min'] = a_mins
    men_data['radius'] = r_men
    men_data['center'] = center
    men_data['zeta'] = zeta
    return men_data


def sinusoidal(target,
               mode='max',
               target_Pc=None,
               surface_tension='pore.surface_tension',
               contact_angle='pore.contact_angle',
               throat_diameter='throat.diameter',
               throat_amplitude='throat.amplitude',
               throat_length='throat.length',
               touch_length='throat.touch_length',
               **kwargs):
    r"""
    The profile of a throat is approximated with a sinusoidal function
    that depends on the average of the connecting pore diameters and throat
    diameter. It represents a converging-diverging geometry that has a minima
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
        'touch' : the maximum capillary pressure a meniscus can sustain before
                  touching a solid feature
        'men' : return the meniscus info for a target pressure
    target_Pc : float (Default is None)
        The target capillary pressure for use with mode 'men'
    surface_tension : dict key (string)
        The dictionary key containing the surface tension values to be used. If
        a pore property is given, it is interpolated to a throat list.
    contact_angle : dict key (string)
        The dictionary key containing the contact angle values to be used. If
        a pore property is given, it is interpolated to a throat list.
    throat_diameter : dict key (string)
        The dictionary key containing the average throat diameter values.
    throat_amplitude : dict key (string)
        The dictionary key containing the amplitude of variation in the throat
        diameter about the mean.
    throat_length : dict key (string)
        The dictionary key containing the throat length values to be used.
    touch_length : dict key (string)
        The dictionary key containing the maximum length that a meniscus can
        protrude into the connecting pore before touching a solid feature and
        therfore invading
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
    # position along throat, amplitude, throat rad, throat length, sigma, theta
    x, A, rt, lt, s, t = syp.symbols('x, A, rt, lt, s, t')
    # Equations
    # Radius profile along throat length
    y = A*syp.cos(2*syp.pi*x) + rt
    # dr/dx used for filling angle
    yprime = y.diff(x)/lt
    # Filling angle
    alpha = syp.atan(yprime)
    # Meniscus Radius of curvature
    R = y/syp.cos(alpha+t)
    # distance from center of curvature to meniscus contact point (Pythagoras)
    a = syp.sqrt(R*R - y*y)
    # angle between throat axis, meniscus center and meniscus contact point
    gamma = syp.atan(y/a)
    # Capillary Pressure function
    f = -2*s*syp.cos(alpha+t)/y
    # Callable expressions
    rx = syp.lambdify((x, A, rt), y, 'numpy')
    Pc = syp.lambdify((x, A, rt, lt, s, t), f, 'numpy')
    rad_curve = syp.lambdify((x, A, rt, lt, s, t), R, 'numpy')
    c2x = syp.lambdify((x, A, rt, lt, s, t), a, 'numpy')
    fill_angle = syp.lambdify((x, A, rt, lt), alpha, 'numpy')
    cap_angle = syp.lambdify((x, A, rt, lt, s, t), gamma, 'numpy')
    theta = np.deg2rad(theta)
    # Network properties
    t_len = network[throat_length]
    pos = np.arange(0.1, 0.9, 1e-3)
    r_amp = network[throat_amplitude]
    r_ts = network[throat_diameter]/2
    # Now find the positions of the menisci along each throat axis
    Y, X = np.meshgrid(r_ts, pos)
    t_Pc = Pc(X, r_amp, Y, t_len, sigma, theta)
    # Values of minima and maxima
    Pc_min = np.min(t_Pc, axis=0)
    Pc_max = np.max(t_Pc, axis=0)
    # Arguments of minima and maxima
    a_min = np.argmin(t_Pc, axis=0)
    a_max = np.argmax(t_Pc, axis=0)
    if mode == 'max':
        return Pc_max
    elif mode == 'touch':
        all_rad = rad_curve(X, r_amp, Y, t_len, sigma, theta)
        all_c2x = c2x(X, r_amp, Y, t_len, sigma, theta)
        all_cen = X*t_len + np.sign(all_rad)*all_c2x
        dist = all_cen + np.abs(all_rad)
        # Only count lengths where meniscus bulges into pore
        dist[all_rad > 0] = 0.0
        touch_len = network[touch_length] + t_len/2
        mask = dist > touch_len
        arg_touch = np.argmax(mask, axis=0)
        # Make sure we only count ones that happen before max pressure
        # And above min pressure (which will be erroneous)
        arg_in_range = (arg_touch < a_max) * (arg_touch > a_min)
        arg_touch[~arg_in_range] = a_max[~arg_in_range]
        x_touch = pos[arg_touch]
        # Return the pressure at which a touch happens
        Pc_touch = Pc(x_touch, r_amp, r_ts, t_len, sigma, theta)
        return Pc_touch
    elif target_Pc is None:
        logger.exception(msg='Please supply a target capillary pressure' +
                         ' when mode is "men"')
    if np.abs(target_Pc) < 1.0:
        target_Pc = 1.0

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
    xpos = pos[arg_x]
    # Output
    men_data = {}
    men_data['pos'] = xpos
    men_data['rx'] = rx(xpos, r_amp, r_ts)
    men_data['alpha'] = fill_angle(xpos, r_amp, r_ts, t_len)
    men_data['c2x'] = c2x(xpos, r_amp, r_ts, t_len, sigma, theta)
    men_data['gamma'] = cap_angle(xpos, r_amp, r_ts, t_len, sigma, theta)
    men_data['radius'] = rad_curve(xpos, r_amp, r_ts, t_len, sigma, theta)
    # xpos is relative to the start of the throat not the center
    # Coop filling assumes center of throat
    men_data['center'] = ((xpos-0.5)*t_len +
                          np.sign(men_data['radius'])*men_data['c2x'])
    logger.info(mode+' calculated for Pc: '+str(target_Pc))
    return men_data
