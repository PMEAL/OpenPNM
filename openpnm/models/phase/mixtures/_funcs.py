import numpy as np
from openpnm.utils import Docorator


docstr = Docorator()


__all__ = [
    'salinity',
    'mixing_rule',
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


def mixing_rule(
    target,
    prop,
    mode='logarithmic',
    power=1,
):
    r"""
    Computes the property of a mixture using the specified mixing rule

    Parameters
    ----------
    target : dict
        The openpnm object to which this model applies
    prop : str
        The dictionary key containing the property of interest on each
        component
    mode : str
        The mixing rule to to use. Options are:

        ============== ========================================================
        mode
        ============== ========================================================
        'logarithmic'  (default) Uses the natural logarithm of the property as:
                       :math:`ln(z) = \Sigma (x_i \cdot ln(\z_i))`
        'linear'       Basic mole fraction weighting of the form
                       :math:`z = \Sigma (x_i \cdot \z_i)`
        'power         Applies an exponent to the property as:
                       :math:`\z^{power} = \Sigma (x_i \cdot \z_i^{power})`
        ============== ========================================================

    power : scalar
        If ``mode='power'`` this indicates the value of the exponent,
        otherwise this is ignored.
    """
    xs = target['pore.mole_fraction']
    ys = target.get_comp_vals(prop)
    z = 0.0
    if mode == 'logarithmic':
        for i in xs.keys():
            z += xs[i]*np.log(ys[i])
        z = np.exp(z)
    elif mode in ['linear', 'simple']:
        for i in xs.keys():
            z += xs[i]*ys[i]
    elif mode == 'power':
        for i in xs.keys():
            z += xs[i]*ys[i]**power
        z = z**(1/power)
    return z


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
