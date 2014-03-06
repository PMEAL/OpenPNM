
"""
module throat_vector
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
    network.set_throat_data(locations=geometry,prop=propname,data=value)

def pore_to_pore(geometry,
                 network,
                 propname,
                 **params):
    r"""
    Calculates throat vector as straight path between connected pores.
    
    Notes
    -----
    There is an important impicit assumption here: the positive direction is
    taken as the direction from the pore with the lower index to the higher.
    This corresponds to the pores in the 1st and 2nd columns of the 
    'connections' array as stored on the network.
    """
    Ts = network.get_throat_indices(geometry.name)
    Ps = network.find_connected_pores(tnums=Ts,flatten=False)
    C0 = network.get_pore_data(prop='coords')[Ps[:,0]]
    C1 = network.get_pore_data(prop='coords')[Ps[:,1]]
    V = C1 - C0
    L = sp.sqrt(sp.sum(V[:,:]**2,axis=1))
    value = V[:,0]/L
    network.set_throat_data(locations=geometry,prop=propname,data=value)
    value = V[:,1]/L
    network.set_throat_data(locations=geometry,prop=propname,data=value)
    value = V[:,2]/L
    network.set_throat_data(locations=geometry,prop=propname,data=value)