import numpy as np
from openpnm.models.phase.mixtures import mixing_rule
from openpnm.models.phase import _phasedocs


__all__ = [
    'gas_mixture_yweighted',
    'gas_pure_TRC',
    'liquid_pure_rp',
    'liquid_mixture_xweighted',
]


@_phasedocs
def liquid_pure_rp(
    phase,
    T='pore.temperature',
    Tc='param.critical_temperature',
    omega='param.acentric_factor',
    Cpg='pore.heat_capacity_gas',
):
    r"""

    Parameters
    ----------
    %(phase)s
    %(T)s
    %(Tc)s
    %(omega)s
    %(Cpg)s

    """
    # Rowlinson and Poling
    T = phase[T]
    Tc = phase[Tc]
    omega = phase[omega]
    Cpgm = phase[Cpg]
    Tr = T/Tc
    if np.any(Tr > 1):
        raise Exception('Cannot calculate liquid property of fluid above'
                        + 'its critical temperature')
    R = 8.314462618
    lhs = 1.586 + 0.49/(1-Tr) \
        + omega*(4.2775 + 6.3*((1-Tr)**(1/3))/Tr + 0.4355/(1-Tr))
    Cp = lhs*R + Cpgm

    return Cp


@_phasedocs
def gas_pure_TRC(
    phase,
    T='pore.temperature',
    a=[],
):
    r"""

    Parameters
    ----------
    %(phase)s
    %(T)s
    a : list
        The coefficients to use (see notes for form of equation). If not
        given the ``phase['param.CAS']`` is used to lookup the values from
        ``chemicals.heat_capacity.TRC_gas_data``

    Returns
    -------

    """
    # TRCCp
    from chemicals.heat_capacity import TRC_gas_data
    T = phase[T]
    if len(a) == 0:
        c = TRC_gas_data.loc[phase.params['CAS']]
        a = list(c[3:11])
    R = 8.314462618
    y = np.zeros_like(T)
    temp = (T - a[7])/(T + a[6])
    mask = T > a[7]
    y[mask] = temp[mask]
    Cp = R*(a[0] + (a[1]/(T**2))*np.exp(-a[1]/T) + a[3]*(y**2)
            + (a[4] - a[5]/((T - a[7])**2))*(y**8))
    return Cp


@_phasedocs
def gas_mixture_yweighted(
    phase,
    Cps='pore.heat_capacity.*',
):
    r"""
    Uses a linearly mole fraction weighted average

    Parameters
    ----------
    %(phase)s
    %(Cps)s

    Returns
    -------

    """
    Cpmix = mixing_rule(phase=phase, prop=Cps, mode='linear')
    return Cpmix


@_phasedocs
def liquid_mixture_xweighted(
    phase,
    Cps='pore.heat_capacity.*',
):
    r"""
    Uses a linearly mole fraction weighted average

    Parameters
    ----------
    %(phase)s
    %(Cps)s

    Returns
    -------

    """
    Cpmix = mixing_rule(phase=phase, prop=Cps, mode='linear')
    return Cpmix
