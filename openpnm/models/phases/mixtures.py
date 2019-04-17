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


def fuller(target,
           molecular_weight='pore.molecular_weight',
           molar_diffusion_volume='pore.molar_diffusion_volume',
           temperature='pore.temperature',
           pressure='pore.pressure'):
    r"""

    """
    A, B = target.components.values()

    T = target[temperature]
    P = target[pressure]
    MA = A[molecular_weight]
    MB = B[molecular_weight]
    vA = A[molar_diffusion_volume]
    vB = B[molar_diffusion_volume]
    MAB = 1e3*2*(1.0/MA + 1.0/MB)**(-1)
    P = P*1e-5
    value = 0.00143*T**1.75/(P*(MAB**0.5)*(vA**(1./3) + vB**(1./3))**2)*1e-4
    return value
