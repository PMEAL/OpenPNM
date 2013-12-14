
"""
module throat_length
===============================================================================

"""
import scipy as sp

def constant(geometry,network,value,**params):
    r"""
    Assigns specified constant value
    """
    network.throat_properties['length'] = value

def straight(geometry,network,**params):
    r"""
    Calculate throat length 
    """
    #Initialize throat_property['length']
    network.throat_properties['length'] = sp.zeros_like(network.throat_properties['type'])
    C1 = network.pore_properties['coords'][network.throat_properties['connections'][:,0]]
    C2 = network.pore_properties['coords'][network.throat_properties['connections'][:,1]]
    E = sp.sqrt(sp.sum((C1-C2)**2,axis=1))  #Euclidean distance between pores
    D1 = network.pore_properties['diameter'][network.throat_properties['connections'][:,0]]
    D2 = network.pore_properties['diameter'][network.throat_properties['connections'][:,1]]
    network.throat_properties['length'] = E - (D1 + D2)/2
        
    