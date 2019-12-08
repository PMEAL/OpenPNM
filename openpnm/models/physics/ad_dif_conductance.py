r"""

.. autofunction:: openpnm.models.physics.ad_dif_conductance.ad_dif

"""

import scipy as _sp


def ad_dif(target,
           conduit_lengths='throat.conduit_lengths',
           pore_pressure='pore.pressure',
           throat_hydraulic_conductance='throat.hydraulic_conductance',
           throat_diffusive_conductance='throat.diffusive_conductance',
           s_scheme='powerlaw'):
    r"""
    Calculate the advective-diffusive conductance of conduits in network, where
    a conduit is ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    conduit_lengths : string
        Dictionary key of the conduit length values

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
        Array containing advective-diffusive conductance values for conduits in
        the geometry attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper conduit_shape_factors argument when
    computig the diffusive and hydraulic conductances.

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
    # Find g for half of pore 1, throat, and half of pore 2
    P = phase[pore_pressure]
    gh = phase[throat_hydraulic_conductance][throats]
    gd = phase[throat_diffusive_conductance][throats]
    gd = _sp.tile(gd, 2)

    Qij = -gh*_sp.diff(P[cn], axis=1).squeeze()
    Qij = _sp.append(Qij, -Qij)

    Peij = Qij / gd
    Peij[(Peij < 1e-10) & (Peij >= 0)] = 1e-10
    Peij[(Peij > -1e-10) & (Peij <= 0)] = -1e-10

    # Export Peclet values (half only since Peij = -Peji)
    phase['throat.peclet.ad'] = _sp.nan
    phase['throat.peclet.ad'][throats] = _sp.absolute(Peij[0:len(Lt)])

    # Correct the flow rate
    Qij = Peij * gd

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
    w = _sp.reshape(w, (throats.size, 2), order='F')
    return w
