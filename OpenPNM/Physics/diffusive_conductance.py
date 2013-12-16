
"""
module diffusive_conductance
===============================================================================

"""

import scipy as sp
import os
propname = os.path.splitext(os.path.basename(__file__))[0]

def constant(physics,network,fluid,value,**params):
    r"""
    Assigns specified constant value
    """
    network.throat_conditions[fluid.name+'_'+propname] = value

def na(physics,network,fluid,**params):
    value = -1
    network.throat_conditions[fluid.name+'_'+propname] = value

def bulk_diffusion(physics,network,fluid,**params):
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
    cp = network.pore_conditions[fluid.name+'_'+'molar_density']
    DABp = network.pore_conditions[fluid.name+'_'+'diffusivity']
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
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    network.throat_conditions[fluid.name+'_'+propname] = value

