import numpy as np
from openpnm.utils import Docorator


docstr = Docorator()

__all__ = [
    "salinity",
    "mole_weighted_average",
    "fuller_diffusivity",
    "wilke_fuller_diffusivity"
]

@docstr.dedent
def salinity(target,
             temperature='pore.temperature',
             concentration='pore.concentration'):
    r"""
    Calculates the salinity in g salt per kg of solution from concentration

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phase.T)s
    concentration : str
        The dictionary key containing the concentration values, in SI units of
        mol/m3.

    Returns
    -------
    salinity : ndarray
        The salinity in g of solute per kg of solution.

    Notes
    -----
    This model is useful for converting known concentration values (e.g.
    calculated by a transport algorithm) into salinity values, which can then
    be used for finding other physical properties of water which are available
    as a function of salinity.

    The salinity correlations are valid for salinity up to 160 g/kg, which
    corresponds to a concentration of 0.05 mol/L (assuming NaCl is the only
    solute)

    """
    C = target[concentration]
    T = target[temperature]
    a = 8.73220929e+00
    b = 6.00389629e+01
    c = -1.19083743e-01
    d = -1.77796042e+00
    e = 3.26987130e-04
    f = -1.09636011e-01
    g = -1.83933426e-07
    S = a + b*C + c*T + d*C**2 + e*T**2 + f*C**3 + g*T**3
    return S


@docstr.dedent
def mole_weighted_average(target, prop):
    r"""

    Parameters
    ----------
    %(models.target.parameters)s
    prop : str
        The dictionary key to the property to be averaged.

    Returns
    -------
    vals : ND-array
        An ND-array containing the mole fraction weighted average value of the
        specified property.
    """
    comps = target.components.values()
    element = prop.split('.')[0]
    if len(comps) == 0:
        vals = np.zeros(target._count(element))*np.nan
    else:
        vals = np.zeros(target._count(element))
        for item in comps:
            frac = target[element + '.mole_fraction.' + item.name]
            temp = item[prop]
            vals += temp*frac
    return vals


def fuller_diffusivity(target, molecular_weight='pore.molecular_weight',
                       molar_diffusion_volume='pore.molar_diffusion_volume',
                       temperature='pore.temperature',
                       pressure='pore.pressure'):
    r"""
    Estimates the diffusion coeffient of both species in a binary gas
    mixture using the Fuller correlation

    Parameters
    ----------
    %(models.target.parameters)s
    molecular_weight : str
        Dictionary key containing the molecular weight of each species.  The
        default is 'pore.molecular_weight'
    molar_diffusion_volume : str
        Dictionary key containing the molar diffusion volume of each species.
        This is used by the Fuller correlation.  The default is
        'pore.molar_diffusion_volume'
    %(models.phase.T)s
    %(models.phase.P)s

    Returns
    -------
    Dij : dict containing ND-arrys
        The dict contains one array for each component, containing the
        diffusion coefficient of that component at each location.
    """
    species_A, species_B = target.components.values()
    T = target[temperature]
    P = target[pressure]
    MA = species_A[molecular_weight]
    MB = species_B[molecular_weight]
    vA = species_A[molar_diffusion_volume]
    vB = species_B[molar_diffusion_volume]
    MAB = 1e3*2*(1.0/MA + 1.0/MB)**(-1)
    P = P*1e-5
    value = 0.00143*T**1.75/(P*(MAB**0.5)*(vA**(1./3) + vB**(1./3))**2)*1e-4
    return value


@docstr.dedent
def wilke_fuller_diffusivity(
        target,
        molecular_weight='pore.molecular_weight',
        molar_diffusion_volume='pore.molar_diffusion_volume',
        temperature='pore.temperature',
        pressure='pore.pressure'):
    r"""
    Estimates the diffusion coeffient of each species in a gas mixture

    Uses the fuller equation to estimate binary diffusivity between pairs, then
    uses the correction of Fairbanks and Wilke to account for the composition
    of the gas mixture.

    Parameters
    ----------
    %(models.target.parameters)s
    molecular_weight : str
        Dictionary key containing the molecular weight of each species.  The
        default is 'pore.molecular_weight'
    molar_diffusion_volume : str
        Dictionary key containing the molar diffusion volume of each species.
        This is used by the Fuller correlation.  The default is
        'pore.molar_diffusion_volume'
    %(models.phase.T)s
    %(models.phase.P)s

    Returns
    -------
    Dij : dict containing ND-arrys
        The dict contains one array for each component, containing the
        diffusion coefficient of that component at each location.

    References
    ----------
    Fairbanks DF and CR Wilke, Diffusion Coefficients in Multicomponent
    Gas Mixtures. Industrial & Engineering Chemistry, 42(3), p471–475 (1950).
    `DOI: 10.1021/ie50483a022 <http://doi.org/10.1021/ie50483a022>`_

    """
    comps = list(target.components.values())
    values = {}
    for i in range(len(comps)):
        A = comps[i]
        denom = 0.0
        for j in range(len(comps)):
            if i != j:
                B = comps[j]
                temp = MixDict(target=target, components=(A, B))
                D = fuller_diffusivity(target=temp,
                                       molecular_weight=molecular_weight,
                                       molar_diffusion_volume=molar_diffusion_volume,
                                       temperature=temperature,
                                       pressure=pressure)
                yB = target['pore.mole_fraction.' + B.name]
                denom += yB/D
        yA = target['pore.mole_fraction.' + A.name]
        values[A.name] = (1 - yA)/denom
    return values


class MixDict(dict):
    r"""
    This utility dict is used to create a temporary mixture object containing
    only two components of a mixture that has several.  This is necessary for
    use of the fuller model for calculating binary diffusivities.
    """
    def __init__(self, target, components):
        super().__init__(target)
        self.components = {}
        for item in components:
            self.components.update({item.name: item})
