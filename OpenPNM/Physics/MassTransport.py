
"""
module MassTransport
===============================================================================

"""

import OpenPNM
import scipy as sp

def DiffusiveConductance(network,fluid):
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
    This function requires that all the necessary transport properties have already been determined.

    """
    cp = fluid['molar_density']
    DABp = fluid['diffusivity']
    ct = network.interpolate_throat_values(cp)
    DABt = network.interpolate_throat_values(DABp)

    gt = ct*DABt*network.throat_properties['diameter']**2/(network.throat_properties['length'])
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    #Find g for half of pore 1
    gp1 = cp*DABp*network.pore_properties['diameter'][pores[:,0]]**2/(network.pore_properties['diameter'][pores[:,0]]/2)
    #Find g for half of pore 2
    gp2 = cp*DABp*network.pore_properties['diameter'][pores[:,1]]**2/(network.pore_properties['diameter'][pores[:,1]]/2)
    g = (1/gt + 1/gp1 + 1/gp2)**(-1)

    return g

