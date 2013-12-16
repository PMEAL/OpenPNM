
"""
module molar_density
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

def ideal_gas(fluid,network,**params):
    r"""
    Uses ideal gas equation to estimate molar density of a pure gas

    """
    R = 8.314
    T = network.pore_conditions[fluid.name+'_'+'temperature']
    P = network.pore_conditions[fluid.name+'_'+'pressure']
    value = P/(R*T)
    network.pore_conditions[fluid.name+'_'+propname] = value