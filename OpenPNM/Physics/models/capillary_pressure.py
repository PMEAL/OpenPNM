r"""
===============================================================================
Submodule -- capillary_pressure
===============================================================================

"""

import scipy as _sp

def butterfly(physics, phase, network, surface_tension='pore.surface_tension',
             contact_angle='pore.contact_angle', throat_diameter='throat.diameter',
             **kwargs):
    r"""
    Computes the capillary entry pressure assuming the throat in a hourglass tube.

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
    The Butterfly equation is:

    .. math::
        P_c = -\frac{2\sigma(cos(arctan(max(dr/dx)) + \theta))}{r}


    """
    fibreRadius = 4.5e-6
    constLength = 1.5e-5
    
    print("Calculating Capillary Pressures...")
    
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
    # Base radius (not including fibre)
    r = network[throat_diameter]/2
    
    r = r[:,_sp.newaxis]
    x = _sp.linspace(0, constLength, 100)
    f = lambda y: fibreRadius*_sp.sin(10*y/(constLength*_sp.pi))
    df = lambda y: (10*fibreRadius/(constLength*_sp.pi))*_sp.cos(10*y/(constLength*_sp.pi))
    r = r - f(x)
    # -2*sigma*cos(theta)/radius
    drdx=_sp.absolute(df(x))
    value = []
    
    for i in range(len(r)):
        if i%1000 == 0:
            print(i)
        radii = r[i]
        caps = []
        deg = theta[i]
        sig = sigma[i]
        caps = []
            
        caps = [-2*sig*_sp.cos(_sp.arctan(drdx[j]) + _sp.radians(deg))/radii[j] for j in range(len(x))]
        maxcap = min(caps)
        
        value.append(maxcap)
    
    '''    
    value = -2*sigma*_sp.cos(_sp.radians(theta))/r
    if throat_diameter.split('.')[0] == 'throat':
        value = value[phase.throats(physics.name)]
    else:
        value = value[phase.pores(physics.name)]
    value[_sp.absolute(value) == _sp.inf] = 0
    '''
    return value

def washburn(physics, phase, network, surface_tension='pore.surface_tension',
             contact_angle='pore.contact_angle', throat_diameter='throat.diameter',
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
    if throat_diameter.split('.')[0] == 'throat':
        value = value[phase.throats(physics.name)]
    else:
        value = value[phase.pores(physics.name)]
    value[_sp.absolute(value) == _sp.inf] = 0
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
            (1 + R/r*(1 - _sp.cos(_sp.radians(alpha)))))
    if throat_diameter.split('.')[0] == 'throat':
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
