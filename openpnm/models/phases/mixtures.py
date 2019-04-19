import numpy as np


def mole_weighted_average(target, prop):
    r"""
    """
    element = prop.split('.')[0]
    vals = np.zeros(target._count(element))
    for item in target.components.keys():
        frac = target[element + '.mole_fraction.' + item]
        temp = target.project[item][prop]
        vals += temp*frac
    return vals


def fuller_diffusivity(target, species_A=None, species_B=None,
                       molecular_weight='pore.molecular_weight',
                       molar_diffusion_volume='pore.molar_diffusion_volume',
                       temperature='pore.temperature',
                       pressure='pore.pressure'):
    r"""

    """
    if (species_A is None) and (species_B is None):
        try:
            species_A, species_B = target.components.values()
        except ValueError:
            raise Exception('fuller_diffusivity applies to binary mixtures only')
    else:
        species_A = target.project[species_A]
        species_B = target.project[species_B]
    T = target[temperature]
    P = target[pressure]
    MA = species_A[molecular_weight]
    MB = species_B[molecular_weight]
    vA = species_A[molar_diffusion_volume]
    vB = species_A[molar_diffusion_volume]
    MAB = 1e3*2*(1.0/MA + 1.0/MB)**(-1)
    P = P*1e-5
    value = 0.00143*T**1.75/(P*(MAB**0.5)*(vA**(1./3) + vB**(1./3))**2)*1e-4
    return value


def wilke_fuller_diffusivity(target, species_A,
                             mole_fraction='pore.mole_fraction',
                             molecular_weight='pore.molecular_weight',
                             molar_diffusion_volume='pore.molar_diffusion_volume',
                             temperature='pore.temperature',
                             pressure='pore.pressure'):
    r"""

    Parameters
    ----------
    target : OpenPNM Mixture object

    Reference
    ---------
    Fairbanks DF and CR Wilke, Diffusion Coefficients in Multicomponent
    Gas Mixtures. Industrial & Engineering Chemistry, 42(3), p471–475 (1950).
    `DOI: 10.1021/ie50483a022 <http://doi.org/10.1021/ie50483a022>`_
    """
    comps = [c for c in target.components.values() if c is not species_A]
    denom = 0.0
    for B in comps:
        D = fuller_diffusivity(target=target,
                               species_A=species_A, species_B=B.name)
        denom += target[mole_fraction + '.' + B.name]/D
    value = (1-target[mole_fraction + '.' + species_A])/(denom)
    return value
