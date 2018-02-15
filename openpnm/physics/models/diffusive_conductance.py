r"""
===============================================================================
Submodule -- diffusive_conductance
===============================================================================

"""

import scipy as _sp
import scipy.constants as _const
import openpnm.utils.misc as misc
import openpnm as pnm


def bulk_diffusion(target, molar_density='pore.molar_density',
                   pore_diffusivity='pore.diffusivity',
                   throat_diffusivity=None,
                   pore_area='pore.area',
                   pore_diameter='pore.diameter',
                   throat_area='throat.area',
                   throat_length='throat.length',
                   throat_diameter='throat.diameter',
                   shape_factor='throat.shape_factor',
                   calc_pore_len=False):
    r"""
    Calculate the diffusive conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas

    Parameters
    ----------
    network : OpenPNM Network Object

    phase : OpenPNM Phase Object
        The phase of interest

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    """
    network = target.simulation.network
    phase = target.simulation.find_phase(target)
    # Get Nt-by-2 list of pores connected to each throat
    Ps = network['throat.conns']
    # Get properties in every pore in the network
    parea = network[pore_area]
    pdia = network[pore_diameter]
    # Get the properties of every throat
    tdia = network[throat_diameter]
    tarea = _sp.pi*(tdia/2)**2
    tlen = network[throat_length]
    # Interpolate pore phase property values to throats
    cp = phase[molar_density]
    ct = phase.interpolate_data(propname=molar_density)
    DABp = phase[pore_diffusivity]
    if throat_diffusivity is None:
        DABt = phase.interpolate_data(propname=pore_diffusivity)
    else:
        DABt = phase[throat_diffusivity]
    if calc_pore_len:
        lengths = misc.conduit_lengths(network, mode='centroid')
        plen1 = lengths[:, 0]
        plen2 = lengths[:, 2]
    else:
        plen1 = (0.5*pdia[Ps[:, 0]])
        plen2 = (0.5*pdia[Ps[:, 1]])
    # Remove any non-positive lengths
    plen1[plen1 <= 0] = 1e-12
    plen2[plen2 <= 0] = 1e-12
    # Find g for half of pore 1
    gp1 = ct*DABt*parea[Ps[:, 0]] / plen1
    gp1[_sp.isnan(gp1)] = _sp.inf
    gp1[~(gp1 > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    # Find g for half of pore 2
    gp2 = ct*DABt*parea[Ps[:, 1]] / plen2
    gp2[_sp.isnan(gp2)] = _sp.inf
    gp2[~(gp2 > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    # Find g for full throat, remove any non-positive lengths
    tlen[tlen <= 0] = 1e-12
    # Get shape factor
    try:
        sf = network[shape_factor]
    except:
        sf = _sp.ones(network.num_throats())
    sf[_sp.isnan(sf)] = 1.0
    gt = (1/sf)*ct*DABt*tarea/tlen
    # Set 0 conductance pores (boundaries) to inf
    gt[~(gt > 0)] = _sp.inf
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[phase.throats(target.name)]
    return value


def mixed_diffusion(target, pore_mixed_diffusivity='pore.mixed_diffusivity',
                    throat_mixed_diffusivity='throat.mixed_diffusivity'):
    r"""
    Uses Knudsen model to adjust the diffusion coefficients to account for the from
    first principles at conditions of interest.

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.
    
    pore_mixed_diffusivity : string
        The dictionary key containing the pore mixed diffusivity values to be used.
    
    throat_mixed_diffusivity : string
        The dictionary key containing the throat mixed diffusivity values to be used.
    
    Notes
    -----
    (1) This model is only valid for dilute systems. Otherwise, the diffusivity
    becomes dependant on concentration and you need to iterate for accuracy.
    
    (2) This model requires `knudsen` model to have already been added to `target`.
    
    """
    phase = target.simulation.find_phase(target)
    # Add `knudsen_scaling` model to phase if `mixed_diffusivity` is not found
    if 'pore.mixed_diffusivity' not in phase.keys():
        knudsen_scaling = pnm.physics.models.diffusion.knudsen_scaling
        target.add_model(propname='pore.mixed_diffusivity',
                         model=knudsen_scaling,
                         diffusivity='pore.diffusivity',
                         knudsen_diffusivity='pore.knudsen_diffusivity')
        target.add_model(propname='throat.mixed_diffusivity',
                         model=knudsen_scaling,
                         diffusivity='throat.diffusivity',
                         knudsen_diffusivity='throat.knudsen_diffusivity')
    target.regenerate_models(propnames=['pore.mixed_diffusivity',
                                        'throat.mixed_diffusivity'])
    # Interleave data to give phase access to `pore/throat.mixed_diffusivity` data
    phase._interleave_data(prop=pore_mixed_diffusivity, sources=[target])
    phase._interleave_data(prop=throat_mixed_diffusivity, sources=[target])
    return bulk_diffusion(target, pore_diffusivity=pore_mixed_diffusivity,
                          throat_diffusivity=throat_mixed_diffusivity)
