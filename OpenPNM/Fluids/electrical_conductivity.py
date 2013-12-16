
"""
module electrical_conductance
===============================================================================

"""
import scipy as sp
import os
propname = os.path.splitext(os.path.basename(__file__))[0]

def constant(fluid,network,value,**params):
    r"""
    Assigns specified constant value
    """
    network.pore_conditions[fluid.name+'_'+propname] = value

def na(fluid,network,**params):
    value = -1
    network.pore_conditions[fluid.name+'_'+propname] = value