r"""
===============================================================================
Submodule -- knudsen
===============================================================================

"""
import scipy.constants as _const
import scipy as _sp


def knudsen(target, diameter='pore.diameter', temperature='pore.temperature', 
            molecular_weight='pore.molecular_weight'):
    r"""
    Computes the pure knudsen diffusivity for.

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.
    
    temperature : string
        The dictionary key containing the temperature values in Kelvin (K)
        
    diameter : string
        The dictionary key containing the pore/throat diameter values to be used.
    
    molecular_weight : float, array_like
        Molecular weight of component A [kg/mol]

    """
    network = target.simulation.network
    phase = target.simulation.find_phase(target)
    Tp = phase[temperature]
    MAp = phase[molecular_weight]
    # Interpolate throat data pores
    Tt = phase.interpolate_data(propname=temperature)
    MAt = phase.interpolate_data(propname=molecular_weight)
    if diameter.split('.')[0] == 'pore':
        DKA = network[diameter]/3 * _sp.sqrt((8*_const.R*Tp)/(_const.pi*MAp))
        DKA = DKA[phase.pores(target.name)]
    elif diameter.split('.')[0] == 'throat':
        DKA = network[diameter]/3 * _sp.sqrt((8*_const.R*Tt)/(_const.pi*MAt))
        DKA = DKA[phase.throats(target.name)]
    else:
        raise Exception('The given diameter is not properly formatted!')
    return DKA


def knudsen_scaling(target, diffusivity='pore.diffusivity',
                    knudsen_diffusivity='pore.knudsen_diffusivity'):
    r"""
    Computes the effective diffusivity considering both bulk and knudsen modes.

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.
    
    diffusivity : string
        The dictionary key containing the pore/throat diffusivity values to be used.
        
    knudsen_diffusivity : string
        The dictionary key containing the pore/throat knudsen diffusivity values
        to be used.
    
    """
    # Verify that the arguments are both pore or throat properties.
    if diffusivity.split('.')[0] != knudsen_diffusivity.split('.')[0]:
        raise Exception('The passed arguments (' + diffusivity + ', ' + \
                        knudsen_diffusivity + ') must be ' + 'both either a ' + \
                        'pore or a throat property!')  
    phase = target.simulation.find_phase(target)
    # Interpolate data for `throat.diffusivity`
    if diffusivity.split('.')[0] == 'throat':
        Db = phase.interpolate_data(propname='pore.'+diffusivity.split('.')[1])
    else:
        Db = phase[diffusivity]
    # Add `knudsen` model to `physics` if not already there.
    if 'pore.knudsen_diffusivity' not in target.keys():
        target.add_model(propname='pore.knudsen_diffusivity', model=knudsen,
                         diameter='pore.diameter')
    if 'throat.knudsen_diffusivity' not in target.keys():
        target.add_model(propname='throat.knudsen_diffusivity', model=knudsen,
                         diameter='throat.diameter')
    target.regenerate_models(propnames=['pore.knudsen_diffusivity',
                                        'throat.knudsen_diffusivity'])
    DK = target[knudsen_diffusivity]
    De = 1 / (1/DK + 1/Db)
    if diffusivity.split('.')[0] == 'pore':
        De = De[phase.pores(target.name)]
    elif diffusivity.split('.')[0] == 'throat':
        De = De[phase.throats(target.name)]
    else:
        raise Exception('The given diameter is not properly formatted!')
    return De
