r"""
===============================================================================
Submodule -- molar_mass
===============================================================================

"""
import scipy as sp

def constant(fluid,network,propname,value,**params):
    r"""
    Assigns specified constant value
    """
    fluid.set_pore_data(prop=propname,data=value)

def na(fluid,network,propname,**params):
    r"""
    Assigns nonsensical, but numerical value of -1.  
    This ensurse stability of other methods 
    but introduces the possibility of being misused.
    """
    value = -1
    fluid.set_pore_data(prop=propname,data=value)

def mixture(fluid,network,propname,MWs=[],mole_fracs=[],**params):
    r"""
    Calculates the average molecular weight of a mixture using mole fraction weighting
    """
    MWs = sp.array(MWs)
    mole_fracs = sp.array(mole_fracs)
    #Ensure mole fraction sum to 1
    fsum = sp.sum(mole_fracs)
    if fsum != 1: print(fluid._logger.warning('mole fractions do not add to 1, so performing normalization'))
    value = sp.sum(MWs*mole_fracs)/fsum
    fluid.set_pore_data(prop=propname,data=value)
