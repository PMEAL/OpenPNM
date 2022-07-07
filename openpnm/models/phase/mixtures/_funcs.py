import numpy as np
from openpnm.utils import Docorator


docstr = Docorator()


__all__ = [
    'salinity',
    'mole_weighted_average',
    'mixture_molecular_weight',
    'mole_summation',
    'from_component',
]


@docstr.dedent
def salinity(
    target,
    temperature='pore.temperature',
    concentration='pore.concentration',
    ):
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


def mixture_molecular_weight(target):
    r"""
    Computes the average molecular weight of a mixture based on mode fractions

    Parameters
    ----------
    %(models.target.parameters)s

    Returns
    -------
    vals : ndarray
        An ND-array containing the mole fraction weighted average molecular
        weight
    """
    xs = [target['pore.mole_fraction.' + c.name] for c in target.components.values()]
    MWs = [c['param.molecular_weight'] for c in target.components.values()]
    MW = np.zeros_like(xs[0])
    for i in range(len(xs)):
        MW += xs[i]*MWs[i]
    return MW


@docstr.dedent
def mole_weighted_average(target, prop):
    r"""
    Computes the mole fraction weighted average of the given property

    Parameters
    ----------
    %(models.target.parameters)s
    prop : string
        The dictionary key to the property to be averaged

    Returns
    -------
    vals : ND-array
        An ND-array containing the mole fraction weighted average value of the
        specified property.
    """
    comps = target.components.values()
    if len(comps) == 0:
        vals = np.zeros(target.Np)*np.nan
    else:
        vals = np.zeros(target.Np)
        for item in comps:
            frac = target['pore.mole_fraction.' + item.name]
            temp = item[prop]
            vals += temp*frac
    return vals


@docstr.dedent
def mole_summation(target):
    r"""
    Computes total mole fraction in each pore given component values

    Parameters
    ----------
    %(models.target.parameters)s

    Returns
    -------
    vals : ND-array
        An ND-array containing the total mole fraction.  Note that this is not
        guaranteed to sum to 1.0.
    """
    xs = [target['pore.mole_fraction.' + c.name] for c in target.components.values()]
    if len(xs) > 0:
        xs = np.sum(xs, axis=0)
    else:
        xs = np.nan
    return xs


@docstr.dedent
def from_component(target, prop, compname):
    r"""
    Fetches the given values from the specified object

    Parameters
    ----------
    %(models.target.parameters)s

    prop : str
        The name of the array to retreive
    compname : str
        The name of the object possessing the desired data

    Returns
    -------
    vals : ND-array
        An ND-array containing the request data
    """
    comp = target.project[compname]
    vals = comp[prop]
    return vals
