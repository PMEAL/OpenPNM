r"""
===============================================================================
Submodule -- molar_mass
===============================================================================

"""
import scipy as _sp

def mixture(fluid,MWs,mole_fracs,**kwargs):
    r"""
    Calculates the average molecular weight of a mixture using mole fraction weighting
    """
    MWs = _sp.array(MWs)
    mole_fracs = _sp.array(mole_fracs)
    #Ensure mole fraction sum to 1
    fsum = _sp.sum(mole_fracs)
    if fsum != 1: print(fluid._logger.warning('mole fractions do not add to 1, so performing normalization'))
    value = _sp.sum(MWs*mole_fracs)/fsum
    return value
