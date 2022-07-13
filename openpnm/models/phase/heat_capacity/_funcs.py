import chemicals as chem
import numpy as np


__all__ = [
    'gas_mixture',
]


def gas_pure(
    target,
    temperature='pore.temperature',
    a=[],
):
    r"""
    """
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
    heat_capacity='pore.heat_capacity',
    mode='linear',
    power=1,
):
    r"""
    """
    Cpmix = target.get_mix_vals(heat_capacity, mode=mode, power=power)
    return Cpmix
