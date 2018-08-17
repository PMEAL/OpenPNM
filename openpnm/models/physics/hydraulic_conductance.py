import scipy as _sp
from scipy import pi


def hagen_poiseuille(target,
                     pore_viscosity='pore.viscosity',
                     throat_viscosity='throat.viscosity',
                     pore_diameter='pore.diameter',
                     throat_diameter='throat.diameter',
                     throat_conduit_lengths='throat.conduit_lengths'):
    r"""
    Calculate the hydraulic conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_viscosity : string
        Dictionary key of the pore viscosity values

    throat_viscosity : string
        Dictionary key of the throat viscosity values

    throat_equivalent_area : string
        Dictionary key of the throat equivalent area values

    throat_conduit_lengths : string
        Dictionary key of the throat conduit lengths

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    # Getting equivalent areas
    D1 = network[pore_diameter][cn[:, 0]]
    Dt = network[throat_diameter][throats]
    D2 = network[pore_diameter][cn[:, 1]]
    # Getting conduit lengths
    L1 = network[throat_conduit_lengths+'.pore1'][throats]
    Lt = network[throat_conduit_lengths+'.throat'][throats]
    L2 = network[throat_conduit_lengths+'.pore2'][throats]
    # Interpolate pore phase property values to throats
    try:
        mut = phase[throat_viscosity]
    except KeyError:
        mut = phase.interpolate_data(propname=pore_viscosity)
    try:
        mup = phase[pore_viscosity]
    except KeyError:
        mup = phase.interpolate_data(propname=throat_viscosity)
    # Find g for half of pore 1
    sf1 = 2*L1/D1
    sf1[sf1 >= 1.0] = 0.999
    gp1 = pi*D1/(64*mup[cn[:, 0]]*_sp.arctanh(sf1))
    gp1[_sp.isnan(gp1)] = _sp.inf
    # Find g for half of pore 2
    sf2 = 2*L2/D2
    sf2[sf2 >= 1.0] = 0.999
    gp2 = pi*D2/(64*mup[cn[:, 1]]*_sp.arctanh(sf2))
    gp2[_sp.isnan(gp2)] = _sp.inf
    # Find g for full throat
    gt = pi*Dt**4/(128*mut[throats]*Lt)
    gt[_sp.isnan(gt)] = _sp.inf
    return (1/gt + 1/gp1 + 1/gp2)**(-1)
