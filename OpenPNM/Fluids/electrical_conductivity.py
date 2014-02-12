
"""
module electrical_conductance
===============================================================================

"""
import scipy as sp

def constant(fluid,network,propname,value,**params):
    r"""
    Assigns specified constant value
    """
    network.set_pore_data(phase=fluid,prop=propname,data=value)

def na(fluid,network,propname,**params):
    value = -1
    network.set_pore_data(phase=fluid,prop=propname,data=value)