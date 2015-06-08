# -*- coding: utf-8 -*-
r"""
===============================================================================
Submodule -- thermal_conductance
===============================================================================

"""
import scipy as sp


def water(phase,
          pore_T='pore.temperature',
          pore_salinity='pore.salinity',
          **kwargs):
    r"""
    Calculates thermal conductivity of pure water or seawater at atmospheric
    pressure using the correlation given in [1]_. Values at temperature higher
    than the normal boiling temperature are calculated at the saturation
    pressure.

    Parameters
    ----------
    phase : OpenPNM Phase Object

    pore_temperature : string
        The dictionary key containing the temperature values.  Temperature must
        be in Kelvin for this emperical equation to work

    pore_salinity : string
        The dictionary key containing the salinity values.  Salinity must be
        expressed in g of salt per kg of solution (ppt).

    Returns
    -------
    k_sw, the thermal conductivity of water/seawater in [W/m.K]

    Notes
    -----
    T must be in K, and S in g of salt per kg of phase, or ppt (parts per
        thousand)
    VALIDITY: 273 < T < 453 K; 0 < S < 160 g/kg;
    ACCURACY: 3 %

    References
    ----------
    [1] D. T. Jamieson, and J. S. Tudhope, Desalination, 8, 393-401, 1970.

    """
    T = phase[pore_T]
    try:
        S = phase[pore_salinity]
    except:
        S = 0
    T68 = 1.00024*T  # convert from T_90 to T_68
    SP = S/1.00472  # convert from S to S_P
    k_sw = 0.001*(10**(sp.log10(240+0.0002*SP) +
                       0.434*(2.3-(343.5+0.037*SP)/T68) *
                       ((1-T68/(647.3+0.03*SP)))**(1/3)))
    value = k_sw
    return value


def chung(phase,
          pore_Cv='pore.heat_capacity',
          pore_acentric='pore.acentric_factor',
          pore_MW='pore.molecular_weight',
          pore_viscosity='pore.viscosity',
          pore_T='pore.temperature',
          pore_Tc='pore.critical_temperature',
          **kwargs):
    r"""
    Uses Chung et al. model to estimate thermal conductivity for gases with
    low pressure(<10 bar) from first principles at conditions of interest

    Parameters
    ----------
    pore_acentric : string
        Dictionary key containing the acentric factor of the component

    pore_Cv :  string
        Dictionary key containing the heat capacity at constant volume
        (J/(mol.K))

    pore_MW : string
        Dictionary key containing the molecular weight of the component
        (kg/mol)

    pore_viscosity : string
        The dictionary key containing the viscosity values (Pa.s)

    pore_T : string
        The dictionary key containing the temperature values (K)

    pore_Tc: string
        The dictionary key containing the critical temperature values (K)

    """
    Cv = phase[pore_Cv]
    acentric = phase[pore_acentric]
    MW = phase[pore_MW]
    R = 8.314
    T = phase[pore_T]
    mu = phase[pore_viscosity]
    Tc = phase[pore_Tc]
    Tr = T/Tc
    z = 2.0 + 10.5*Tr**2
    beta = 0.7862 - 0.7109*acentric + 1.3168*acentric**2
    alpha = Cv/R - 3/2
    s = 1 + alpha*((0.215+0.28288*alpha-1.061*beta+0.26665*z) /
                   (0.6366+beta*z+1.061*alpha*beta))
    value = 3.75*s*(mu)*R/(MW)
    return value


def sato(phase,
         pore_MW='pore.molecular_weight',
         pore_Tb='pore.boiling_point',
         pore_T='pore.temperature',
         pore_Tc='pore.critical_temperature',
         **params):
    r"""
    Uses Sato et al. model to estimate thermal conductivity for pure liquids
    from first principles at conditions of interest

    Parameters
    ----------
    pore_Tb :  string
        Dictionary key containing the toiling temperature of the component (K)

    pore_MW : string
        Dictionary key containing the molecular weight of the component
        (kg/mol)

    pore_T : string
        The dictionary key containing the temperature values (K)

    pore_Tc: string
        The dictionary key containing the critical temperature values (K)

    """
    T = phase[pore_T]
    Tc = phase[pore_Tc]
    MW = phase[pore_MW]
    Tbr = phase[pore_Tb]/Tc
    Tr = T/Tc
    value = (1.11/((MW*1e3)**0.5))*(3+20*(1-Tr)**(2/3))/(3+20*(1-Tbr)**(2/3))
    return value
