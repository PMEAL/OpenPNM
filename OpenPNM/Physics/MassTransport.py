
"""
module MassTransport
===============================================================================

"""

import scipy as sp

def DiffusiveConductance(network,fluid):
    r"""
    Calculate the diffusive conductance of conduits in network ( 1/2 pore - full throat - 1/2 pore ) based on the area

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
        The fluid of interest

    Notes
    -----
    This function requires that all the necessary fluid properties should already be calculated.

    """
    try:
        cp = fluid.pore_conditions['molar_density']
        DABp = fluid.pore_conditions['diffusivity']
    except:
        raise Exception('Necessary fluid properies are not present ' + fluid['name'])
    ct = fluid.interpolate_throat_conditions(network,cp)
    DABt = fluid.interpolate_throat_conditions(network,DABp)

    #Get Nt-by-2 list of pores connected to each throat
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    #Find g for half of pore 1
    gp1 = ct*DABt*network.pore_properties['diameter'][pores[:,0]]**2/(0.5*network.pore_properties['diameter'][pores[:,0]])
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = ct*DABt*network.pore_properties['diameter'][pores[:,1]]**2/(0.5*network.pore_properties['diameter'][pores[:,1]])
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    gt = ct*DABt*network.throat_properties['diameter']**2/(network.throat_properties['length'])
    g = (1/gt + 1/gp1 + 1/gp2)**(-1)
    fluid.throat_conditions['diffusive_conductance'] = g

