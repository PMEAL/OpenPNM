
"""
module vapor_pressure
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

def Antoine(fluid,network,A=8.07131,B=1730.63,C=233.426,**params):
    r"""
    Uses Antoine equation to estimate vapor pressure of a pure component

    Parameters
    ----------
    A, B, C :  float, array_like
            Antoine vapor pressure constants for pure compounds

    """
    T = network.get_pore_conditions(fluid.name,'temperature')
    value = (10**(A-B/(C+T-273.15)))*1.333e2
    network.set_pore_condition(fluid.name,propname,value)