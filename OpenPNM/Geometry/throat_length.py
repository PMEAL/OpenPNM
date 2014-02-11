
"""
module throat_length
===============================================================================

"""
import scipy as sp
import os
propname = os.path.splitext(os.path.basename(__file__))[0]
propname = propname.split('_')[1]

def constant(geometry,network,value,**params):
    r"""
    Assigns specified constant value
    """
    network.set_throat_data(labels=geometry,prop=propname,data=value)

def straight(geometry,network,**params):
    r"""
    Calculate throat length 
    """
    #Initialize throat_property['length']
    network.set_throat_data(prop='length',data=sp.zeros_like(network.get_throat_indices(indices=False)))
    C1 = network.get_pore_data(prop='coords')[network.get_throat_data(prop='connections')[:,0]]
    C2 = network.get_pore_data(prop='coords')[network.get_throat_data(prop='connections')[:,1]]
    E = sp.sqrt(sp.sum((C1-C2)**2,axis=1))  #Euclidean distance between pores
    D1 = network.get_pore_data(prop='diameter')[network.get_throat_data(prop='connections')[:,0]]
    D2 = network.get_pore_data(prop='diameter')[network.get_throat_data(prop='connections')[:,1]]
    value = E-(D1+D2)/2
    network.set_throat_data(labels=geometry,prop=propname,data=value)
        
    