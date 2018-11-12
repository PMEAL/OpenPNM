from .misc import generic_conductance


def hagen_poiseuille(target,
                     pore_area='pore.area',
                     throat_area='throat.area',
                     pore_viscosity='pore.viscosity',
                     throat_viscosity='throat.viscosity',
                     conduit_lengths='throat.conduit_lengths',
                     conduit_shape_factors='throat.flow_shape_factors'):
    r"""
    Calculate the hydraulic conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

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

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    conduit_shape_factors : string
        Dictionary key of the conduit FLOW shape factor values

    Returns
    -------
    g : ndarray
        Array containing hydraulic conductance values for conduits in the
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
    return generic_conductance(target=target, transport_type='flow',
                               pore_area=pore_area,
                               throat_area=throat_area,
                               pore_diffusivity=pore_viscosity,
                               throat_diffusivity=throat_viscosity,
                               conduit_lengths=conduit_lengths,
                               conduit_shape_factors=conduit_shape_factors)


def hagen_poiseuille_2D(target,
                        pore_diameter='pore.diameter',
                        throat_diameter='throat.diameter',
                        pore_viscosity='pore.viscosity',
                        throat_viscosity='throat.viscosity',
                        conduit_lengths='throat.conduit_lengths',
                        conduit_shape_factors='throat.flow_shape_factors'):
    r"""
    Returns hydraulic conductance for a 2D pore-throat-pore assembly (conduit).

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    # Getting pore/throat diameters
    D1 = network[pore_diameter][cn[:, 0]]
    Dt = network[throat_diameter][throats]
    D2 = network[pore_diameter][cn[:, 1]]
    # Getting conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors+'.pore1'][throats]
        SFt = phase[conduit_shape_factors+'.throat'][throats]
        SF2 = phase[conduit_shape_factors+'.pore2'][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Getting viscosity values
    mut = phase[throat_viscosity][throats]
    mu1 = phase[pore_viscosity][cn[:, 0]]
    mu2 = phase[pore_viscosity][cn[:, 1]]
    # Find g for half of pore 1, throat, and half of pore 2
    g1 = D1**3 / (12*mu1*L1)
    g2 = D2**3 / (12*mu2*L2)
    gt = Dt**3 / (12*mut*Lt)

    return (1/gt/SFt + 1/g1/SF1 + 1/g2/SF2)**(-1)
