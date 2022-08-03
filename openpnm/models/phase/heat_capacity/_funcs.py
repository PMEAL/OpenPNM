import numpy as np
from openpnm.models.phase.mixtures import mixing_rule


__all__ = [
    'gas_mixture_yweighted',
    'gas_pure_TRC',
    'liquid_pure_rp',
    'liquid_mixture_xweighted',
]


def liquid_pure_rp(
    target,
    T='pore.temperature',
    Tc='param.critical_temperature',
    omega='param.acentric_factor',
    Cpg='pore.heat_capacity_gas',
):
    r"""
    """
    # Rowlinson and Poling
    T = target[T]
    Tc = target[Tc]
    omega = target[omega]
    Cpgm = target[Cpg]
    Tr = T/Tc
    if np.any(Tr > 1):
        raise Exception('Cannot calculate liquid property of fluid above'
                        + 'its critical temperature')
    R = 8.314462618
    lhs = 1.586 + 0.49/(1-Tr) \
        + omega*(4.2775 + 6.3*((1-Tr)**(1/3))/Tr + 0.4355/(1-Tr))
    Cp = lhs*R + Cpgm

    return Cp


def gas_pure_TRC(
    target,
    T='pore.temperature',
    a=[],
):
    r"""
    """
    # TRCCp
    from chemicals.heat_capacity import TRC_gas_data
    T = target[T]
    if len(a) == 0:
        c = TRC_gas_data.loc[target.params['CAS']]
        a = list(c[3:11])
    R = 8.314462618
    y = np.zeros_like(T)
    temp = (T - a[7])/(T + a[6])
    mask = T > a[7]
    y[mask] = temp[mask]
    Cp = R*(a[0] + (a[1]/(T**2))*np.exp(-a[1]/T) + a[3]*(y**2)
            + (a[4] - a[5]/((T - a[7])**2))*(y**8))
    return Cp


def gas_mixture_yweighted(
    target,
    Cps='pore.heat_capacity.*',
):
    r"""
    """
    Cpmix = mixing_rule(target=target, prop=Cps, mode='linear')
    return Cpmix


def liquid_mixture_xweighted(
    target,
    Cps='pore.heat_capacity.*',
):
    r"""
    """
    Cpmix = mixing_rule(target=target, prop=Cps, mode='linear')
    return Cpmix
