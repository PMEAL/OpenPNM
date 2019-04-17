r"""

.. autofunction:: openpnm.models.physics.diffusive_conductance.dispersion
.. autofunction:: openpnm.models.physics.diffusive_conductance.generic_conductance

"""

import scipy as _sp


def dispersion(target,
               pore_area='pore.area',
               throat_area='throat.area',
               pore_diffusivity='pore.diffusivity',
               throat_diffusivity='throat.diffusivity',
               conduit_lengths='throat.conduit_lengths',
               conduit_shape_factors='throat.poisson_shape_factors',
               pore_pressure='pore.pressure',
               throat_hydraulic_conductance='throat.hydraulic_conductance',
               throat_diffusive_conductance='throat.diffusive_conductance',
               s_scheme='powerlaw'):
    r"""
    Calculate the advective-diffusive conductance of conduits in network,
    where a conduit is ( 1/2 pore - full throat - 1/2 pore ).
    See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    pore_diffusivity : string
        Dictionary key of the pore diffusivity values

    throat_diffusivity : string
        Dictionary key of the throat diffusivity values

    conduit_lengths : string
        Dictionary key of the conduit length values

    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values

    pore_pressure : string
        Dictionary key of the pore pressure values

   throat_hydraulic_conductance : string
       Dictionary key of the throat hydraulic conductance values

   throat_diffusive_conductance : string
       Dictionary key of the throat diffusive conductance values

   s_scheme : string
       Name of the space discretization scheme to use

    Returns
    -------
    g : ndarray
        Array containing dispersive conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper flow_shape_factor argument.

    """
    return generic_conductance(
        target=target,
        transport_type='dispersion',
        pore_area=pore_area,
        throat_area=throat_area,
        pore_diffusivity=pore_diffusivity,
        throat_diffusivity=throat_diffusivity,
        conduit_lengths=conduit_lengths,
        conduit_shape_factors=conduit_shape_factors,
        pore_pressure=pore_pressure,
        throat_hydraulic_conductance=throat_hydraulic_conductance,
        throat_diffusive_conductance=throat_diffusive_conductance,
        s_scheme=s_scheme)


def generic_conductance(target, transport_type, pore_area, throat_area,
                        pore_diffusivity, throat_diffusivity,
                        conduit_lengths, conduit_shape_factors, **kwargs):
    r"""
    Calculate the generic conductance (could be mass, thermal, electrical,
    ionic, or hydraylic) of conduits in the network, where a conduit is
    ( 1/2 pore - full throat - 1/2 pore ).

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    transport_type : string
        Dictionary key of the transport type

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    pore_diffusivity : string
        Dictionary key of the pore diffusivity values

    throat_diffusivity : string
        Dictionary key of the throat diffusivity values

    conduit_lengths : string
        Dictionary key of the conduit length values

    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values

    Returns
    -------
    g : ndarray
        Array containing conductance values for conduits in the geometry
        attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper shape factor.

    (4) shape_factor depends on the physics of the problem, i.e. diffusion-like
    processes and fluid flow need different shape factors.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    # Getting conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    # Preallocating g
    g1, g2, gt = _sp.zeros((3, len(Lt)))
    # Setting g to inf when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    g1[~m1] = g2[~m2] = gt[~mt] = _sp.inf
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors+'.pore1'][throats]
        SFt = phase[conduit_shape_factors+'.throat'][throats]
        SF2 = phase[conduit_shape_factors+'.pore2'][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Find g for half of pore 1, throat, and half of pore 2
    if transport_type == 'dispersion':
        for k, v in kwargs.items():
            if k == 'pore_pressure':
                pore_pressure = v
            elif k == 'throat_hydraulic_conductance':
                throat_hydraulic_conductance = v
            elif k == 'throat_diffusive_conductance':
                throat_diffusive_conductance = v
            elif k == 's_scheme':
                s_scheme = v

        P = phase[pore_pressure]
        gh = phase[throat_hydraulic_conductance]
        gd = phase[throat_diffusive_conductance]
        gd = _sp.tile(gd, 2)

        Qij = -gh*_sp.diff(P[cn], axis=1).squeeze()
        Qij = _sp.append(Qij, -Qij)

        Peij = Qij/gd
        Peij[(Peij < 1e-10) & (Peij >= 0)] = 1e-10
        Peij[(Peij > -1e-10) & (Peij <= 0)] = -1e-10
        Qij = Peij*gd

        if s_scheme == 'upwind':
            w = gd + _sp.maximum(0, -Qij)
        elif s_scheme == 'hybrid':
            w = _sp.maximum(0, _sp.maximum(-Qij, gd-Qij/2))
        elif s_scheme == 'powerlaw':
            w = gd * _sp.maximum(0, (1 - 0.1*_sp.absolute(Peij))**5) + \
                _sp.maximum(0, -Qij)
        elif s_scheme == 'exponential':
            w = -Qij / (1 - _sp.exp(Peij))
        else:
            raise Exception('Unrecognized discretization scheme: ' + s_scheme)
        w = _sp.reshape(w, (network.Nt, 2), order='F')
        return w
    else:
        raise Exception('Unknown keyword for "transport_type", can only be' +
                        ' "dispersion"')
    # Apply shape factors and calculate the final conductance
    return (1/gt/SFt + 1/g1/SF1 + 1/g2/SF2)**(-1)
