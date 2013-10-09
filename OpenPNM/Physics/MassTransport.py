
"""
module MassTransport
===============================================================================

"""

import OpenPNM
import scipy as sp


def Binary_Diff_In_a_Half_Pore_Throat_Half_Pore_Conduit(network,c=0,DAB=0):
    r"""
    This does not account for temperature gradients
    """
    gt = c*DAB*network.throat_properties['diameter']**2/network.throat_properties['length']
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    gp1 = c*DAB*network.pore_properties['diameter'][pores[:,0]]**2/(network.pore_properties['diameter'][pores[:,0]]/2)
    gp2 = c*DAB*network.pore_properties['diameter'][pores[:,1]]**2/(network.pore_properties['diameter'][pores[:,1]]/2) 
    g = (1/gt + 1/gp1 + 1/gp2)**(-1)
    return g
    
def Binary_Diff_In_a_Half_Pore_Throat_Half_Pore_Conduit_with_T_Dependence(network,T,P):
    r"""
    T is a vector, or T already exists in pore_properties
    """
    c = P/R/network.pore_properties['T']
    DAB = fuller_equation(T,P)
    gt = c*DAB*network.throat_properties['diameter']**2/network.throat_properties['length']
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    gp1 = c*DAB*network.pore_properties['diameter'][pores[:,0]]**2/(network.pore_properties['diameter'][pores[:,0]]/2)
    gp2 = c*DAB*network.pore_properties['diameter'][pores[:,1]]**2/(network.pore_properties['diameter'][pores[:,1]]/2) 
    g = (1/gt + 1/gp1 + 1/gp2)**(-1)
    return g