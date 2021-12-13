from chemicals import numba_vectorized
import chemicals as chem


def liquid_heat_capacity(
    target,
    Cpgm='pore.gas_heat_capacity',
    T='pore.temperature'
):
    r"""
    Computes the heat capacity of a liquid using the ``Rowlinson_Poling``
    method from ``chemicals.heat_capacity``.

    Parameters
    ----------

    """
    T = target['pore.temperature']
    Cpgm = target['pore.gas_heat_capacity']
    Tc = target['param.critical_temperature']
    omega = target['param.acentric_factor']
    Cplm = numba_vectorized.Rowlinson_Poling(T, Tc, omega, Cpgm)
    return Cplm


def gas_heat_capacity(
    target,
    temperature='pore.temperature'
):
    r"""
    Computes the heat capacity of a gas using the ``TRCCp`` method
    from ``chemicals.heat_capacity``.

    Parameters
    ----------

    """
    T = target[temperature]
    props = chem.heat_capacity.TRC_gas_data.loc[target.settings['CAS']]
    _, Tmin, Tmax, a0, a1, a2, a3, a4, a5, a6, a7, I, J, Hfg = props
    Cp = numba_vectorized.TRCCp(T, a0, a1, a2, a3, a4, a5, a6, a7)
    return Cp