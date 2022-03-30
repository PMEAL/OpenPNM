from openpnm.utils import Docorator


docstr = Docorator()

__all__ = ["standard", "ideal_gas", "water"]

@docstr.get_sections(base='models.phase.density', sections=['Returns'])
@docstr.dedent
def standard(target,
             mol_weight='pore.molecular_weight',
             molar_density='pore.molar_density'):
    r"""
    Calculates the mass density from the molecular weight and molar density

    Parameters
    ----------
    %(models.target.parameters)s
    mol_weight : str
        Name of the dictionary key on ``target`` where the array containing
        molecular weight values is stored
    molar_density : str
        Name of the dictionary key on ``target`` where the array containing
        molar density values is stored

    Returns
    -------
    densities : ndarray
        A numpy ndarray containing density values

    Notes
    -----
    %(models.phase.SI_note)s

    """
    MW = target[mol_weight]
    rho = target[molar_density]
    value = rho*MW
    return value


@docstr.dedent
def ideal_gas(target, pressure='pore.pressure', temperature='pore.temperature',
              mol_weight='pore.molecular_weight'):
    r"""
    Uses ideal gas law to calculate the mass density of an ideal gas

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phase.T)s
    %(models.phase.P)s
    mol_weight : str
        Name of the dictionary key on ``target`` where the array containing
        molecular weight values is stored

    Returns
    -------
    %(models.phase.density.returns)s

    """

    P = target[pressure]
    T = target[temperature]
    MW = target[mol_weight]
    R = 8.31447
    value = P/(R*T)*MW
    return value


def water(target, temperature='pore.temperature', salinity='pore.salinity'):
    r"""
    Calculates density of pure water or seawater at atmospheric pressure
    using Eq. (8) given by Sharqawy et. al [1]. Values at temperature higher
    than the normal boiling temperature are calculated at the saturation
    pressure.

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phase.T)s
    salinity : str
        Name of the dictionary key on ``target`` where the array containing
        salinity values is stored.

    Returns
    -------
    %(models.phase.density.returns)s

    Notes
    -----
     T must be in K, and S in g of salt per kg of phase, or ppt (parts per
        thousand)
    VALIDITY: 273 < T < 453 K; 0 < S < 160 g/kg;
    ACCURACY: 0.1 %

    References
    ----------
    [1] Sharqawy M. H., Lienhard J. H., and Zubair, S. M., Desalination and
    Water Treatment, 2010.

    """
    T = target[temperature]
    if salinity in target.keys():
        S = target[salinity]
    else:
        S = 0
    a1 = 9.9992293295E+02
    a2 = 2.0341179217E-02
    a3 = -6.1624591598E-03
    a4 = 2.2614664708E-05
    a5 = -4.6570659168E-08
    b1 = 8.0200240891E-01
    b2 = -2.0005183488E-03
    b3 = 1.6771024982E-05
    b4 = -3.0600536746E-08
    b5 = -1.6132224742E-11
    TC = T-273.15
    rho_w = a1 + a2*TC + a3*TC**2 + a4*TC**3 + a5*TC**4
    d_rho = b1*S + b2*S*TC + b3*S*(TC**2) + b4*S*(TC**3) + b5*(S**2)*(TC**2)
    rho_sw = rho_w + d_rho
    value = rho_sw
    return value
