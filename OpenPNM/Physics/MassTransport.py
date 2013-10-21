
"""
module MassTransport
===============================================================================

"""

import OpenPNM
import scipy as sp

def DiffusiveConductance(network,fluid_name):
    r"""
    Calculate the diffusive conductance of conduits in network ( 1/2 pore - full throat - 1/2 pore ) based on the area

    Parameters
    ----------
    network : OpenPNM Network Object
        The network for which conductance should be calculated

    fluid : OpenPNM Fluid Object
        The fluid of interest

    Notes
    -----
    This function requires that all the necessary fluid properties have already been determined.

    """
    cp = network.pore_conditions['molar_density'+'_'+fluid_name]
    DABp = network.pore_conditions['diffusivity'+'_'+fluid_name]
    ct = network.interpolate_throat_values(cp)
    DABt = network.interpolate_throat_values(DABp)

    #Get Nt-by-2 list of pores connected to each throat
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    #Find g for half of pore 1
    gp1 = cp*DABp*network.pore_properties['diameter'][pores[:,0]]**2/(network.pore_properties['diameter'][pores[:,0]]/2)
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = cp*DABp*network.pore_properties['diameter'][pores[:,1]]**2/(network.pore_properties['diameter'][pores[:,1]]/2)
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    gt = ct*DABt*network.throat_properties['diameter']**2/(network.throat_properties['length'])
    g = (1/gt + 1/gp1 + 1/gp2)**(-1)
    network.throat_conditions['diffusive_conductance'+'_'+fluid_name] = g

