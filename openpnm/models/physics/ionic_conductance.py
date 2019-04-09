import scipy as _sp


def poisson(target,
            pore_area='pore.area',
            throat_area='throat.area',
            conduit_lengths='throat.conduit_lengths',
            conduit_shape_factors='throat.poisson_shape_factors'):
    r"""
    """
    return _generic_conductance(target=target,
                                transport_type='poisson',
                                pore_area=pore_area,
                                throat_area=throat_area,
                                pore_diffusivity='',
                                throat_diffusivity='',
                                conduit_lengths=conduit_lengths,
                                pore_volume='',
                                conduit_shape_factors=conduit_shape_factors)


def laplace(target,
            pore_area='pore.area',
            throat_area='throat.area',
            pore_diffusivity='pore.diffusivity',
            throat_diffusivity='throat.diffusivity',
            conduit_lengths='throat.conduit_lengths',
            conduit_shape_factors='throat.poisson_shape_factors',
            pore_volume='pore.volume',
            pore_temperature='pore.temperature',
            throat_temperature='throat.temperature',
            pore_valence='pore.valence',
            throat_valence='throat.valence',
            pore_concentration='pore.concentration'):
    r"""
    """
    return _generic_conductance(target=target,
                                transport_type='laplace',
                                pore_area=pore_area,
                                throat_area=throat_area,
                                pore_diffusivity=pore_diffusivity,
                                throat_diffusivity=throat_diffusivity,
                                conduit_lengths=conduit_lengths,
                                pore_volume=pore_volume,
                                pore_temperature=pore_temperature,
                                throat_temperature=throat_temperature,
                                pore_valence=pore_valence,
                                throat_valence=throat_valence,
                                pore_concentration=pore_concentration,
                                conduit_shape_factors=conduit_shape_factors)


def electroneutrality(target,
                      pore_area='pore.area',
                      throat_area='throat.area',
                      pore_diffusivity='pore.diffusivity',
                      throat_diffusivity='throat.diffusivity',
                      conduit_lengths='throat.conduit_lengths',
                      conduit_shape_factors='throat.poisson_shape_factors',
                      pore_volume='pore.volume',
                      pore_temperature='pore.temperature',
                      throat_temperature='throat.temperature',
                      pore_valence='pore.valence',
                      throat_valence='throat.valence',
                      pore_concentration='pore.concentration'):
    r"""
    """
    return _generic_conductance(target=target,
                                transport_type='electroneutrality',
                                pore_area=pore_area,
                                throat_area=throat_area,
                                pore_diffusivity=pore_diffusivity,
                                throat_diffusivity=throat_diffusivity,
                                conduit_lengths=conduit_lengths,
                                pore_volume=pore_volume,
                                pore_temperature=pore_temperature,
                                throat_temperature=throat_temperature,
                                pore_valence=pore_valence,
                                throat_valence=throat_valence,
                                pore_concentration=pore_concentration,
                                conduit_shape_factors=conduit_shape_factors)


def _generic_conductance(target, transport_type, pore_area, throat_area,
                         pore_diffusivity, throat_diffusivity, conduit_lengths,
                         pore_volume, pore_temperature, throat_temperature,
                         pore_valence, throat_valence, pore_concentration,
                         conduit_shape_factors, **kwargs):
    r"""
    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    # Getting equivalent areas
    A1 = network[pore_area][cn[:, 0]]
    At = network[throat_area][throats]
    A2 = network[pore_area][cn[:, 1]]
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
    # Poisson or Laplace
    if transport_type in ['poisson', 'laplace']:
        g1[m1] = (A1)[m1] * L1[m1]
        g2[m2] = (A2)[m2] * L2[m2]
        gt[mt] = (At)[mt] * Lt[mt]
    # Electroneutrality
    elif transport_type == 'electroneutrality':
        F = 96485.3329
        R = 8.3145
        # Getting pores volumes
        Vol1 = network[pore_volume][cn[:, 0]]
        Vol2 = network[pore_volume][cn[:, 1]]
        # Interpolate pore phase property values to throats
        try:
            Tt = phase[throat_temperature][throats]
        except KeyError:
            Tt = phase.interpolate_data(propname=pore_temperature)[throats]
        try:
            T1 = phase[pore_temperature][cn[:, 0]]
            T2 = phase[pore_temperature][cn[:, 1]]
        except KeyError:
            T1 = phase.interpolate_data(propname=throat_temperature)[cn[:, 0]]
            T2 = phase.interpolate_data(propname=throat_temperature)[cn[:, 1]]
        # Finding species present in the ionic solution
        ions = []
        for key in phase.keys():
            if key[:len(pore_diffusivity)] == pore_diffusivity:
                ions.append(key[len(pore_diffusivity):])
        for i in ions:
            # Check if a concetration field is defined
            try:
                c1 = phase[pore_concentration+i][cn[:, 0]]
                c2 = phase[pore_concentration+i][cn[:, 1]]
            except KeyError:
                c1 = _sp.zeros((network.Nt))[cn[:, 0]]
                c2 = _sp.zeros((network.Nt))[cn[:, 1]]
            ct = (c1*Vol1 + c2*Vol2)/(Vol1 + Vol2)
        # Interpolate pore phase property values to throats
            try:
                Dt = phase[throat_diffusivity+i][throats]
                Vt = phase[throat_valence+i][throats]
            except KeyError:
                Dt = phase.interpolate_data(
                    propname=pore_diffusivity+i)[throats]
                Vt = phase.interpolate_data(
                    propname=pore_valence+i)[throats]
            try:
                D1 = phase[pore_diffusivity+i][cn[:, 0]]
                D2 = phase[pore_diffusivity+i][cn[:, 1]]
                V1 = phase[pore_valence+i][cn[:, 0]]
                V2 = phase[pore_valence+i][cn[:, 1]]
            except KeyError:
                D1 = phase.interpolate_data(
                    propname=throat_diffusivity+i)[cn[:, 0]]
                D2 = phase.interpolate_data(
                    propname=throat_diffusivity+i)[cn[:, 1]]
                V1 = phase.interpolate_data(
                    propname=throat_valence+i)[cn[:, 0]]
                V2 = phase.interpolate_data(
                    propname=throat_valence+i)[cn[:, 1]]

            g1[m1] += F**2 * V1**2 * (D1*A1*c1)[m1] / (R * T1 * L1[m1])
            g2[m2] += F**2 * V2**2 * (D2*A2*c2)[m1] / (R * T2 * L2[m2])
            gt[mt] += F**2 * Vt**2 * (Dt*At*ct)[mt] / (R * Tt * Lt[mt])

    else:
        raise Exception('Unknown keyword for "transport_type", can only be' +
                        ' "poisson", "laplace" or "electroneutrality"')
    # Preallocating g_inv
    g_inv1, g_inv2, g_invt = _sp.zeros((3, len(Lt)))
    f1, f2, ft = [gi != 0 for gi in [g1, g2, gt]]
    g_inv1[f1] = 1/g1[f1]
    g_inv2[f2] = 1/g2[f2]
    g_invt[ft] = 1/gt[ft]
    # Apply shape factors and calculate the final conductance
    return (g_inv1/SF1 + g_inv2/SF2 + g_invt/SFt)**(-1)
