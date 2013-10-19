
"""
module MassTransport
===============================================================================

"""

import OpenPNM
import scipy as sp

def DiffusiveConductance(network):
    r"""
    Calculate the diffusive conductance of conduits in network ( 1/2 pore - full throat - 1/2 pore )

    Parameters
    ----------
    network : OpenPNM Network Object

    Notes
    -----
    This function requires that all the necessary transport properties have already been determined.

    """
    R = 8.314462
    P = network.pore_conditions['pressure']
    T = network.pore_conditions['temperature']
    network.pore_conditions['molar_density'] = P/(R*T)
    network.interpolate_throat_conditions('molar_density')
    ct = network.throat_conditions['molar_density']
    cp = network.pore_conditions['molar_density']
    network.interpolate_throat_conditions('diffusion coefficient')
    DABt = network.throat_conditions['diffusion coefficient']
    DABp = network.pore_conditions['diffusion coefficient']

    gt = ct*DABt*network.throat_properties['diameter']**2/network.throat_properties['length']
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    #Find g for half of pore 1
    gp1 = cp*DABp*network.pore_properties['diameter'][pores[:,0]]**2/(network.pore_properties['diameter'][pores[:,0]]/2)
    #Find g for half of pore 2
    gp2 = cp*DABp*network.pore_properties['diameter'][pores[:,1]]**2/(network.pore_properties['diameter'][pores[:,1]]/2)
    g = (1/gt + 1/gp1 + 1/gp2)**(-1)

    return g

