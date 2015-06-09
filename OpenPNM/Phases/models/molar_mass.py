r"""
===============================================================================
Submodule -- molar_mass
===============================================================================

"""
import scipy as sp


def mixture(phase,
            molar_mass='pore.molar_mass',
            mole_frac='pore.mole_fraction',
            **kwargs):
    r"""
    Calculates the average molecular weight of a mixture using mole fraction
    weighting

    Parameters
    ----------
    phase : OpenPNM Phase Object
        The phase for which the molar mass is to be calculated.  This phase
        must have sub-phases.

    molar_mass : string, optional (default = 'pore.molar_mass')
        The property name for the molar masses of each sub-phase

    mole_frac : string, optional (default = 'pore.mole_fraction')
        The property name for the mole fraction of each sub-phase
    """
    MW = sp.zeros((phase.Np,))
    for comp in phase._phases:
        MW = MW + comp[molar_mass]*comp[mole_frac]
    return MW
