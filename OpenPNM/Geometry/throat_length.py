r"""
===============================================================================
Submodule -- throat_length
===============================================================================

"""
import scipy as sp

def constant(geometry,
             network,
             propname,
             value,
             **params):
    r"""
    Assigns specified constant value
    """
    network.set_data(prop=propname,throats=geometry.throats(),data=value)

def straight(geometry,
             network,
             propname,
             pore_diameter='diameter',
             **params):
    r"""
    Calculate throat length 
    """
    #Initialize throat_property['length']
    pore1 = network.get_data(prop='conns',throats=geometry.throats())[:,0]
    pore2 = network.get_data(prop='conns',throats=geometry.throats())[:,1]
    C1 = network.get_data(prop='coords',pores=pore1)
    C2 = network.get_data(prop='coords',pores=pore2)
    E = sp.sqrt(sp.sum((C1-C2)**2,axis=1))  #Euclidean distance between pores
    D1 = network.get_data(prop=pore_diameter,pores=pore1)
    D2 = network.get_data(prop=pore_diameter,pores=pore2)
    value = E-(D1+D2)/2
    network.set_data(prop=propname,throats=geometry.throats(),data=value)
        
    