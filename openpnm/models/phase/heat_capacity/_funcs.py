import chemicals as chem
import numpy as np
from openpnm.models.phase.mixtures import mixing_rule


__all__ = [
    'gas_mixture',
    'gas_pure',
    'liquid_pure',
    'liquid_mixture',
]


def liquid_pure(
    target,
    temperature='pore.temperature',
    critical_temperature='param.critical_temperature',
    acentric_factor='param.acentric_factor',
    heat_capacity_gas='pore.heat_capacity_gas',
):
    r"""
    """
    # Rowlinson and Poling
    T = target[temperature]
    Tc = target[critical_temperature]
    omega = target[acentric_factor]
    Cpgm = target[heat_capacity_gas]
    Tr = T/Tc
    if np.any(Tr > 1):
        raise Exception('Cannot calculate liquid property of fluid above'
                        + 'its critical temperature')
    R = 8.314462618
    lhs = 1.586 + 0.49/(1-Tr) \
        + omega*(4.2775 + 6.3*((1-Tr)**(1/3))/Tr + 0.4355/(1-Tr))
    Cp = lhs*R + Cpgm

    return Cp


def gas_pure(
    target,
    temperature='pore.temperature',
    a=[],
):
    r"""
    """
    # TRCCp
    T = target[temperature]
    if len(a) == 0:
        c = chem.heat_capacity.TRC_gas_data.loc[target.params['CAS']]
        a = list(c[3:11])
    R = 8.314462618
    y = np.zeros_like(T)
    temp = (T - a[7])/(T + a[6])
    mask = T > a[7]
    y[mask] = temp[mask]
    Cp = R*(a[0] + (a[1]/(T**2))*np.exp(-a[1]/T) + a[3]*(y**2)
            + (a[4] - a[5]/((T - a[7])**2))*(y**8))
    return Cp


def gas_mixture(
    target,
    heat_capacity='pore.heat_capacity.*',
):
    r"""
    """
    Cpmix = mixing_rule(target=target, prop=heat_capacity, mode='linear')
    return Cpmix


def liquid_mixture(
    target,
    heat_capacity='pore.heat_capacity.*',
):
    r"""
    """
    Cpmix = mixing_rule(target=target, prop=heat_capacity, mode='linear')
    return Cpmix
