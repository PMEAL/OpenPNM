r"""
===============================================================================
Submodule -- molar_mass
===============================================================================

"""
import scipy as sp

def mixture(phase,MWs,mole_fracs,**kwargs):
    r"""
    Calculates the average molecular weight of a mixture using mole fraction weighting
    """
    MWs = sp.array(MWs)
    mole_fracs = sp.array(mole_fracs)
    #Ensure mole fraction sum to 1
    fsum = sp.sum(mole_fracs)
    if fsum != 1: print(phase._logger.warning('mole fractions do not add to 1, so performing normalization'))
    value = sp.sum(MWs*mole_fracs)/fsum
    return value
