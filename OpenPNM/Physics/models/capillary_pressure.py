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
            (1 + R/r*(1 - _sp.cos(_sp.radians(alpha)))))
    value = value[phase.throats(physics.name)]
    return value

def static_pressure(network,
                    physics,
                    phase,
                    density = 'pore.density',
                    throat_prop = 'throat.occupancy',
                    g = [0,0,9.81],
                    **kwargs):
    r'''
    Finds the highest point on each cluster and adds the corresponding static
    fluid pressure to the entry pressure of each throat.

    Parameters
    ----------
    throat_prop : dictionary key (string)
        The name of the array on the phase object describing the phase
        distribution.  This can be a bit tricky to calculate, see examples
        for more details.

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
    It is important to remember that the 'top' of the Network corresponds to
    the maximum coordinate value.  The static pressure is thus calculated using
    the distance from the 'top' of the Network, and not coordiante value
    directly.

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.Cubic(shape=[25,1,25], spacing=0.0001)
    >>> water = OpenPNM.Phases.Water(network=pn)
    >>> water['pore.density'] = 997  # kg/m3
    >>> phys_water = OpenPNM.Physics.GenericPhysics(network=pn,
                                                    phase=water,
                                                    pores=pn.Ps,
                                                    throats=pn.Ts)

    Add the 'static_pressure' model to the water Physics object:

    >>> phys_water.models.add(model=OpenPNM.Physics.models.capillary_pressure.static_pressure,
                              propname='pore.static_pressure',
                              throat_prop='throat.occupancy',
                              density='pore.density',
                              regen_mode='deferred')

    Rigorously, it is necessary to create an IP algorithm to determine a water
    distribution in the Network, but for the sake of this example, an
    artificial distribution will be used:

    >>> water['pore.occupancy'] = sp.rand(pn.Np,) < 0.3

    It is necessary to find the throats that define clusters on the phase of
    interest:

    >>> Ts = pn.find_neighbor_throats(water['pore.occupancy'], mode='intersection')
    >>> water['throat.occupancy'] = False
    >>> water['throat.occupancy'][Ts] = True
    >>> phys_water.models.regenerate()

    To visualize the result use:

    .. code-block:: python

        plt.matshow(pn.asarray(phys_water['pore.static_pressure'])[:,0,:].T,
                    interpolation='none',
                    origin='lower')

    '''
    # Setup model variables and parameters
    static_pressure = _sp.zeros((network.Np,))
    rho = phase[density]
    g = _sp.array(g)
    # Labels clusters of defending phase
    clusters = network.find_clusters(physics[throat_prop])
    # Set non-invaded pores to label of -1
    Ps = network.find_connected_pores(physics[throat_prop], flatten=True)
    Ps = network.tomask(Ps)
    clusters[~Ps] = -1
    # Remove the -1 cluster from list
    cluster_nums = _sp.unique(clusters)
    cluster_nums = cluster_nums[~_sp.in1d(cluster_nums,-1)]
    # Scan through each labelled cluster and find static pressure within
    for cluster in cluster_nums:
        Ps = _sp.where(clusters == cluster)[0]
        tops = _sp.amax(network['pore.coords'][Ps,:], axis=0)
        h = tops - network['pore.coords'][Ps]
        P_temp = g*h
        P_temp = _sp.reshape(P_temp[:, _sp.where(g>0)[0]], -1)
        static_pressure[Ps] = P_temp*rho[Ps]
    return static_pressure

























