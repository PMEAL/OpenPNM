
"""
module hydraulic_conductance
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

def hagen_poiseuille(physics,network,fluid,**params):
    r"""
    Calculates the hydraulic conductvity of throat assuming square geometry using a modified Hagen Poiseuille model

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
    """
    mup = network.pore_conditions[fluid.name+'_'+'viscosity']
    mut = fluid.interpolate_throat_conditions(network,mup)
    #Get Nt-by-2 list of pores connected to each throat
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    #Find g for half of pore 1
    gp1 = 2.28*(network.pore_properties['diameter'][pores[:,0]]/2)**4/(network.pore_properties['diameter'][pores[:,0]]*mut)
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = 2.28*(network.pore_properties['diameter'][pores[:,1]]/2)**4/(network.pore_properties['diameter'][pores[:,1]]*mut)
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    gt = 2.28*(network.throat_properties['diameter']/2)**4/(2*network.throat_properties['length']*mut)
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    network.throat_conditions[fluid.name+'_'+propname] = value



