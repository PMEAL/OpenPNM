
"""
module MassTransport
===============================================================================

"""

import OpenPNM
import scipy as sp

def Chapman_Enskog_DiffusionCoefficient(T,P):


    return

def Fuller(T,P,MA=31.99,MB=28.01,VA=16.6,VB=17.9):

    MAB = 2*(1/MA+1/MB)**(-1)
    D = 0.00143*T**1.75/(P*MAB**0.5*(VA**(1/3)+VB**(1/3))**2)*1e-4

    return D

def FullerScaling(network,DABo=2.09e-5,To=298.,Po=101325.,Ti=298.,Pi=101325.):
    r"""
    Uses the Fuller model to adjust a diffusion coeffient from reference conditions to conditions of interest

    Parameters
    ----------
    network : OpenPNM Network Object

    DABo : float, array_like
        Diffusion coefficient at reference conditions

    Po, To, Pi & Ti : float, array_like
        Pressure & temperature at reference conditions, pressure & temperature at conditions of interest, respectively

    """
    try:
        Ti = network.pore_conditions['temperature']
    except:
        network.pore_conditions['temperature'] = Ti
    try:
        Pi = network.pore_conditions['pressure']
    except:
        network.pore_conditions['pressure'] = Pi
    Di = DABo*(Po/To**1.75)*(Ti**1.75/Pi)
    network.pore_conditions['diffusion coefficient'] = Di
    #Return network for GUI
    return {'network': network}

def Conduits_DiffusionConductivity(network):
    r"""
    Calculate the effective diffusive conductivity of a 1/2 pore - throat - 1/2 pore conduit

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

