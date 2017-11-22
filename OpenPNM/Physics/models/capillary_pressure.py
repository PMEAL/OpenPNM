r"""
===============================================================================
Submodule -- capillary_pressure
===============================================================================

"""

import scipy as _sp
import numpy as np
from OpenPNM.Base import logging
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
    if (surface_tension.split('.')[0] == 'pore' and
       diameter.split('.')[0] == 'throat'):
        sigma = phase[surface_tension]
        sigma = phase.interpolate_data(data=sigma)
    else:
        sigma = phase[surface_tension]
    if (contact_angle.split('.')[0] == 'pore' and
       diameter.split('.')[0] == 'throat'):
        theta = phase[contact_angle]
        theta = phase.interpolate_data(data=theta)
    else:
        theta = phase[contact_angle]
    return element, sigma, theta


def _handle_zeros(array, mode='max', value=None):
    r"""
    Convert zeros in an array to either the max, min or specified value
    Useful for handling pores or throats with zero diameter i.e. boundaries

    Parameters
    ----------
    mode : Determines what value to replace zeros with, uses non-zero values.
    options are max, min , mean
    """
    if value is None:
        if mode == 'max':
            value = array.max()
        elif mode == 'min':
            value = array[~array == 0.0].min()
        elif mode == 'mean':
            value = array[~array == 0.0].mean()
    array[array == 0.0] = value
    return array


def washburn(physics, phase, network, surface_tension='pore.surface_tension',
             contact_angle='pore.contact_angle', diameter='throat.diameter',
             **kwargs):
    r"""
    Computes the capillary entry pressure assuming the throat in a cylindrical tube.

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
    diameter : dict key (string)
        The dictionary key containing the element diameter values to be used.

    Notes
    -----
    The Washburn equation is:

    .. math::
        P_c = -\frac{2\sigma(cos(\theta))}{r}

    This is the most basic approach to calculating entry pressure and is
    suitable for highly non-wetting invading phases in most materials.

    """
    element, sigma, theta = _get_key_props(phase=phase,
                                           diameter=diameter,
                                           surface_tension=surface_tension,
                                           contact_angle=contact_angle)
    r = network[diameter]/2
    # Take care of any zeros - Boundary pores should be invaded with ease
    r = _handle_zeros(r, mode='max')
    value = -2*sigma*_sp.cos(_sp.radians(theta))/r
    if element == 'throat':
        value = value[phase.throats(physics.name)]
    else:
        value = value[phase.pores(physics.name)]
    value[_sp.absolute(value) == _sp.inf] = 0
    return value


def purcell(physics, phase, network, r_toroid,
            surface_tension='pore.surface_tension',
            contact_angle='pore.contact_angle',
            diameter='throat.diameter',
            **kwargs):
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
    diameter : dict key (string)
        The dictionary key containing the element diameter values to be used.
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

    element, sigma, theta = _get_key_props(phase=phase,
                                           diameter=diameter,
                                           surface_tension=surface_tension,
                                           contact_angle=contact_angle)
    r = network[diameter]/2
    # Take care of any zeros - Boundary pores should be invaded with ease
    r = _handle_zeros(r, mode='max')
    R = r_toroid
    alpha = theta - 180 + _sp.arcsin(_sp.sin(_sp.radians(theta)/(1+r/R)))
    value = (-2*sigma/r) * \
        (_sp.cos(_sp.radians(theta - alpha)) /
            (1 + R/r*(1 - _sp.cos(_sp.radians(alpha)))))
    if element == 'throat':
        value = value[phase.throats(physics.name)]
    else:
        value = value[phase.pores(physics.name)]
    return value


def purcell_bi(physics, phase, network, r_toroid,
            surface_tension='pore.surface_tension',
            contact_angle='pore.contact_angle',
            diameter='throat.diameter',
            h_max='pore.diameter',
            max_dist=True,
            **kwargs):
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
    entity = diameter.split('.')[0]
    if surface_tension.split('.')[0] == 'pore' and entity == 'throat':
        sigma = phase[surface_tension]
        sigma = phase.interpolate_data(data=sigma)
    else:
        sigma = phase[surface_tension]
    if contact_angle.split('.')[0] == 'pore' and entity == 'throat':
        theta = phase[contact_angle]
        theta = phase.interpolate_data(data=theta)
    else:
        theta = phase[contact_angle]
    theta = 180 - theta #Mason and Morrow have the definitions switched
    th = _sp.deg2rad(theta)
    rt = network[diameter]/2
    R = r_toroid
    a_min = th -np.pi + np.arcsin((np.sin(th))/(1+rt/R))
    a_max = th - np.arcsin((np.sin(th))/(1+rt/R))
    if max_dist and entity == 'throat':
        #Perform analysis for entry into both pores
        r_max = np.zeros([network.Nt,2])
        for j in range(2):
            Pj = network['throat.conns'][:,j]
            dj = network[h_max][Pj]
            max_reached = np.zeros(network.Nt, dtype=bool)
            alpha_reached = np.zeros(network.Nt)
            "With increasing filling angle assess whether interface has passed the"
            "critical distance and record the critical angle"
            a_space = _sp.linspace(1e-3, _sp.pi, 181)
            nudge = 0.001
            for a_test in a_space:
                nudgers = network.throats()[np.around(th-a_test,0) == 90]
                if len(nudgers) > 0:
                    th[nudgers]+=nudge
                r = R*(1+(rt/R)-_sp.cos(a_test))/_sp.cos(th-a_test)
                #Vertical adjustment for centre of circle
                y_off = R*np.sin(a_test)
                #Angle between contact point - centre - vertical
                zeta = (th-a_test-(np.pi/2))
                c = y_off - r*np.cos(zeta)
                y_max = c+r
                ts = network.throats()[(y_max > (dj)) * (~max_reached)]
                if len(ts) > 0:
                    max_reached[ts]=True
                    alpha_reached[ts]=a_test
                if len(nudgers) > 0:
                    th[nudgers]-=nudge
            "Any interfaces that never reach a wall are ok"
            alpha_reached[~max_reached] = a_max[~max_reached]
            print("Percentage max before BT",100*np.sum(max_reached[_sp.absolute(a_max)>_sp.absolute(alpha_reached)])/len(a_max))
            "Any interfaces that can expand to maximum curvature before hitting a"
            "pore wall are ok"
            alpha_reached[_sp.absolute(a_max)<_sp.absolute(alpha_reached)] = \
                          a_max[_sp.absolute(a_max)<_sp.absolute(alpha_reached)]
            r_max[:,j] = R*(1+(rt/R)-_sp.cos(alpha_reached))/_sp.cos(th-alpha_reached)
            #r_max = R*(1+(rt/R)-_sp.cos(alpha_reached))/_sp.cos(th-alpha_reached)

        value= (2*np.vstack((sigma,sigma)).T)/r_max
    else:
        r_max = R*(1+(rt/R)-_sp.cos(a_max))/_sp.cos(th-a_max)
        value = 2*sigma/r_max

    if entity == 'throat':
        value = value[phase.throats(physics.name)]
    else:
        value = value[phase.pores(physics.name)]
    return value

def static_pressure(network,
                    physics,
                    phase,
                    pore_density='pore.density',
                    pore_occupancy='pore.occupancy',
                    g=[0, 0, 9.81],
                    **kwargs):
    r'''
    Finds the highest point on each cluster and adds the corresponding static
    fluid pressure to the entry pressure of each throat.

    Parameters
    ----------
    pore_occupancy : dictionary key (string)
        The name of the array on the phase object describing the phase
        distribution.

    density : dictionary key (string)
        String providing the dictionary location of the phase density.  The
        default is 'pore.density'.

    g : list
        A three component vector describing the direction and magnitude of the
        force acting on the fluid.  The default is [0,0,9.81] corresponding to
        Earth's gravity acting in the downward z-direction.

    Returns
    -------
    An Np long list containing the static fluid pressure within each pore.

    Notes
    -----
    (1) It is important to remember that the 'top' of the Network corresponds
    to the maximum coordinate value.  The static pressure is thus calculated
    using the distance from the 'top' of the Network.

    (2) There is a slight flaw in the logic of sending the pore occupancy,
    rather than throat occupancy: cluster labeling using pore occupancy invokes
    site percolation rather then bond percolation.  Hence, although it is
    physically possible for two neighboring pores to be on different clusters,
    this method will count them as on the same clusters.  This inaccuracy was
    necessary, however, so that the method worked for both defending and
    invading phase.

    Examples
    --------
    >>> import OpenPNM
    >>> import scipy as sp
    >>> pn = OpenPNM.Network.Cubic(shape=[25,1,50], spacing=0.0001)
    >>> water = OpenPNM.Phases.Water(network=pn)
    >>> water['pore.density'] = 997  # kg/m3
    >>> phys_water = OpenPNM.Physics.GenericPhysics(network=pn,
    ...                                             phase=water,
    ...                                             pores=pn.Ps,
    ...                                             throats=pn.Ts)

    Add the 'static_pressure' model to the water Physics object:

    >>> f = OpenPNM.Physics.models.capillary_pressure.static_pressure
    >>> phys_water.models.add(model=f,
    ...                       propname='pore.static_pressure',
    ...                       pore_occupancy='pore.occupancy',
    ...                       density='pore.density',
    ...                       regen_mode='deferred')

    Rigorously speaking, it is necessary to create an IP algorithm to determine
    a water distribution in the Network, but for the sake of this example, an
    artificial distribution will be used:

    >>> water['pore.occupancy'] = sp.rand(pn.Np,) < 0.5
    >>> phys_water.models.regenerate()

    To visualize the result use:

    .. code-block:: python

        plt.matshow(pn.asarray(phys_water['pore.static_pressure'])[:,0,:].T,
                    interpolation='none',
                    origin='lower')

    '''
    # Setup model variables and parameters
    static_pressure = _sp.zeros((network.Np,))
    rho = phase[pore_density]
    g = _sp.array(g)
    # Labels clusters of defending phase
    clusters = network.find_clusters2(phase[pore_occupancy])
    # Remove the -1 cluster from list
    cluster_nums = _sp.unique(clusters)
    cluster_nums = cluster_nums[~_sp.in1d(cluster_nums, -1)]
    # Scan through each labelled cluster and find static pressure within
    for cluster in cluster_nums:
        Ps = _sp.where(clusters == cluster)[0]
        tops = _sp.amax(network['pore.coords'][Ps, :], axis=0)
        h = tops - network['pore.coords'][Ps]
        P_temp = g*h
        P_temp = _sp.reshape(P_temp[:, _sp.where(g > 0)[0]], -1)
        static_pressure[Ps] = P_temp*rho[Ps]
    return static_pressure


def cuboid(physics, phase, network,
           surface_tension='pore.surface_tension',
           contact_angle='pore.contact_angle',
           diameter='throat.diameter', **kwargs):
    r"""
    Computes the capillary entry pressure assuming the throat in a cube tube.

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
    diameter : dict key (string)
        The dictionary key containing the element diameter values to be used.

    Notes
    -----
    The equation is taken from Non-equilibrium effects in capillarity and
    interfacial area in two-phase flow: dynamic pore-network modelling

    """
    element, sigma, theta = _get_key_props(phase=phase,
                                           diameter=diameter,
                                           surface_tension=surface_tension,
                                           contact_angle=contact_angle)
    # Convert theta to rad
    theta *= 2*_sp.pi/360
    rad = network[diameter]/2
    # Take care of any zeros - Boundary pores should be invaded with ease
    rad = _handle_zeros(rad, mode='max')
    Theta = ((theta+_sp.cos(theta)**2-_sp.pi/4-_sp.sin(theta)*_sp.cos(theta)) /
             (_sp.cos(theta)-_sp.sqrt(_sp.pi/4-theta+_sp.sin(theta) *
              _sp.cos(theta))))
    value = (sigma/rad)*Theta
    if element == 'throat':
        value = value[phase.throats(physics.name)]
    else:
        value = value[phase.pores(physics.name)]
    value[_sp.absolute(value) == _sp.inf] = 0
    return value


def from_throat(physics, phase, network,
                capillary_pressure='throat.capillary_pressure',
                operator='min',
                **kwargs):
    r"""
    The capillary pressure for a pore is calculated from the adjoining throats

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network object is
    phase : OpenPNM Phase Object
        Phase object for the invading phases containing the surface tension and
        contact angle values.
    capillary_pressure : string
        label for throat data to use
    operator : string
        Operator for throat values to convert to pore value
        Accepted values are min, max, mean
    """
    value = np.zeros(network.Np)
    functions = {'min': np.min,
                 'max': np.max,
                 'mean': np.mean}
    if operator not in functions.keys():
        operator = 'mean'

    for i in range(network.Np):
        ts = network.find_neighbor_throats(pores=i)
        value[i] = functions[operator](physics[capillary_pressure][ts])

    return value


def kelvin(physics, phase, network, diameter='pore.diameter',
           temperature='pore.temperature',
           vapor_pressure='pore.vapor_pressure',
           molecular_weight='pore.molecular_weight',
           density='pore.density',
           surface_tension='pore.surface_tension',
           **kwargs):
    r"""
    Calculate the critical vapor pressure that causes droplets to condense or
    evaporate inside a pore. Only works with site percolation
    """

    T = phase[temperature]
    P0 = phase[vapor_pressure]
    M = phase[molecular_weight]
    rho = phase[density]
    gamma = phase[surface_tension]
    r = network[diameter]/2
    # Take care of any zeros - Boundary pores should be invaded with ease
    r = _handle_zeros(r, mode='max')
    R = 8.314
    value = P0*np.exp((M*2*gamma)/(rho*R*T*r))
    return value

def ransohoff_snap_off(physics, phase, network,
                       shape_factor=2.0,
                       require_pair=True,
                       contact_angle='pore.contact_angle',
                       surface_tension='pore.surface_tension',
                       throat_diameter='throat.diameter',
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
    network : OpenPNM Network Object
        The Network object is
    phase : OpenPNM Phase Object
        Phase object for the invading phases containing the surface tension and
        contact angle values.
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
    if (surface_tension.split('.')[0] == 'pore'):
        sigma = phase[surface_tension]
        sigma = phase.interpolate_data(data=sigma)
    else:
        sigma = phase[surface_tension]
    if (contact_angle.split('.')[0] == 'pore'):
        theta = phase[contact_angle]
        theta = phase.interpolate_data(data=theta)
    else:
        theta = phase[contact_angle]
    try:
        geometry = network.geometries(network.geometries()[0])[0]
        all_verts = geometry[vertices]
    except:
        logger.error("Model will only work if geometry has property: "+
                     vertices)
    # Work out whether throat geometry can support at least one pair of
    # adjacent arc menisci that can grow and merge to form snap-off
    # Only works if throat vertices are in convex hull order
    angles_ok = np.zeros(network.Nt, dtype=bool)
    for T in range(network.Nt):
        verts = all_verts[T]
        x = verts[:,0]
        y = verts[:,1]
        z = verts[:,2]
        verts_plus = np.vstack((np.roll(x,1), np.roll(y,1), np.roll(z,1))).T
        verts_minus = np.vstack((np.roll(x,-1), np.roll(y,-1), np.roll(z,-1))).T
        v1 = verts_plus - verts
        v2 = verts_minus - verts
        corner_angles = np.rad2deg(tr.angle_between_vectors(v1, v2, axis=1))
        # Logical test for existence of arc menisci
        am = theta[T] <= 90 - corner_angles/2
        if require_pair:
            # Logical test for two adjacent arc menisci
            am_pair = np.any(np.logical_or(np.logical_and(am, np.roll(am, +1)),
                                           np.logical_and(am, np.roll(am, -1))))
            angles_ok[T] = am_pair
        else:
            # Logical test for any arc menisci
            angles_ok[T] = np.any(am)
        
    # Condition for arc menisci to form in corners 
    rad_Ts = network[throat_diameter]/2
    # Ransohoff and Radke eq. 4
    C = 1/rad_Ts - 1/wavelength
    value = sigma*C
    # Only throats that can support arc menisci can snap-off
    value[~angles_ok] = np.nan
    logger.info("Snap off pressures calculated for "
                +str(np.around(100*np.sum(angles_ok)/np.size(angles_ok),0))
                + "% of throats")
    return value[phase.throats(physics.name)]

