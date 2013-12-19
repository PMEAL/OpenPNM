
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
    network.set_pore_condition(fluid.name,propname,value)

def na(fluid,network,**params):
    value = -1
    network.set_pore_condition(fluid.name,propname,value)