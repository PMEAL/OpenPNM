
"""
module throat_length
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
    propname = propname.split('_')[1] #remove leading pore_ or throat_ from dictionary key
    network.set_throat_data(locations=geometry,prop=propname,data=value)

def straight(geometry,
             network,
             propname,
             pore_diameter='diameter',
             **params):
    r"""
    Calculate throat length 
    """
    #Initialize throat_property['length']
    C1 = network.get_pore_data(prop='coords')[network.get_throat_data(prop='connections')[:,0]]
    C2 = network.get_pore_data(prop='coords')[network.get_throat_data(prop='connections')[:,1]]
    E = sp.sqrt(sp.sum((C1-C2)**2,axis=1))  #Euclidean distance between pores
    D1 = network.get_pore_data(prop=pore_diameter)[network.get_throat_data(prop='connections')[:,0]]
    D2 = network.get_pore_data(prop=pore_diameter)[network.get_throat_data(prop='connections')[:,1]]
    value = E-(D1+D2)/2
    network.set_throat_data(locations=geometry,prop=propname,data=value)
        
    