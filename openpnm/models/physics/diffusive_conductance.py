import openpnm as op
import scipy as _sp


def ordinary_diffusion(target,
                       pore_diffusivity='pore.diffusivity',
                       throat_diffusivity='throat.diffusivity',
                       throat_equivalent_area='throat.equivalent_area',
                       throat_conduit_lengths='throat.conduit_lengths'):
    r"""
    Calculate the diffusive conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diffusivity : string
        Dictionary key of the pore diffusivity values

    throat_diffusivity : string
        Dictionary key of the throat diffusivity values

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

    """
    network = target.project.network
    phase = target.project.find_phase(target)
    geom = target.project.find_geometry(target)
    cn = network['throat.conns']
    # Getting equivalent areas
    A1 = geom[throat_equivalent_area+'.pore1']      # Equivalent area pore 1
    At = geom[throat_equivalent_area+'.throat']     # Equivalent area throat
    A2 = geom[throat_equivalent_area+'.pore2']      # Equivalent area pore 2
    # Getting conduit lengths
    L1 = geom[throat_conduit_lengths+'.pore1']       # Equivalent length pore 1
    Lt = geom[throat_conduit_lengths+'.throat']      # Equivalent length throat
    L2 = geom[throat_conduit_lengths+'.pore2']       # Equivalent length pore 2
    # Interpolate pore phase property values to throats
    try:
        DABt = phase[throat_diffusivity]
    except KeyError:
        DABt = phase.interpolate_data(propname=pore_diffusivity)
    try:
        DABp = phase[pore_diffusivity]
    except KeyError:
        DABp = phase.interpolate_data(propname=throat_diffusivity)
    # Remove any non-positive lengths
    L1[L1 <= 0] = 1e-12
    L2[L2 <= 0] = 1e-12
    Lt[Lt <= 0] = 1e-12
    # Find g for half of pore 1
    gp1 = DABp[cn[:, 0]]*A1 / L1
    gp1[_sp.isnan(gp1)] = _sp.inf
    gp1[gp1<=0] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    # Find g for half of pore 2
    gp2 = DABp[cn[:, 1]]*A2 / L2
    gp2[_sp.isnan(gp2)] = _sp.inf
    gp2[gp2<=0] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    # Find g for full throat
    gt = DABt*At / Lt
    gt[gt<=0] = _sp.inf
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[phase.throats(target.name)]
    return value
