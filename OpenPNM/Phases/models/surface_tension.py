r"""
===============================================================================
Submodule -- surface_tension
===============================================================================

"""

import scipy as sp


def water(phase,
          pore_T='pore.temperature',
          pore_salinity='pore.salinity',
          **kwargs):
    r"""
    Calculates surface tension of pure water or seawater at atmospheric
    pressure using Eq. (28) given by Sharqawy et. al [1]_. Values at
    temperature higher than the normal boiling temperature are calculated at
    the saturation pressure.

    Parameters
    ----------
    pore_temperature : string
        The dictionary key containing the temperature values.  Temperature must
        be in Kelvin for this emperical equation to work

    pore_salinity : string
        The dictionary key containing the salinity values.  Salinity must be
        expressed in g of salt per kg of solution (ppt).

    Returns
    -------
    sigma_sw, the surface tension of seawater in [N/m]

    Notes
    -----
     T must be in K, and S in g of salt per kg of phase, or ppt (parts per
        thousand)
    VALIDITY: 273 < T < 313 K; 0 < S < 40 g/kg;
    ACCURACY: 0.2 %

    References
    ----------
    [1] Sharqawy M. H., Lienhard J. H., and Zubair, S. M., Desalination and
    Water Treatment, 2010.

    """
    T = phase['pore.temperature']
    try:
        S = phase['pore.salinity']
    except:
        S = 0
    sigma_w = 0.2358*((1-(T/647.096))**1.256)*(1-0.625*(1-(T/647.096)))
    a1 = 2.2637334337E-04
    a2 = 9.4579521377E-03
    a3 = 3.3104954843E-02
    TC = T-273.15
    sigma_sw = sigma_w*(1+(a1*TC+a2)*sp.log(1+a3*S))
    value = sigma_sw
    return value


def eotvos(phase, k,
           pore_T='pore.temperature',
           pore_Tc='pore.critical_temperature',
           pore_molar_density='pore.molar_density', **kwargs):
    r"""
    Missing description

    Parameters
    ----------
    k : float
        Constant parameter specific to fluid

    pore_T : string
        The dictionary key containing the temperature values (K)

    pore_Tc : string
        The dictionary key containing the critical temperature values (K)

    pore_molar_density : string
        The dictionary key containing the molar density values (K)


    TODO: Needs description, and improve definition of k

    """
    Tc = phase[pore_Tc]
    T = phase[pore_T]
    Vm = 1/phase[pore_molar_density]
    value = k*(Tc-T)/(Vm**(2/3))
    return value


def guggenheim_katayama(phase, K2, n,
                        pore_T='pore.temperature',
                        pore_Tc='pore.critical_temperature',
                        pore_Pc='pore.critical_pressure',
                        **kwargs):
    r"""
    Missing description

    Parameters
    ----------
    K2 : scalar
        Fluid specific constant

    n : scalar
        Fluid specific constant

    pore_T : string
        The dictionary key containing the temperature values (K)

    pore_Tc : string
        The dictionary key containing the critical temperature values (K)

    pore_Pc : string
        The dictionary key containing the critical pressure values (K)

    TODO: Needs description
    """
    T = phase[pore_T]
    Pc = phase[pore_Pc]
    Tc = phase[pore_Tc]
    sigma_o = K2*Tc**(1/3)*Pc**(2/3)
    value = sigma_o*(1-T/Tc)**n
    return value


def brock_bird_scaling(phase, sigma_o, To,
                       pore_T='pore.temperature',
                       pore_Tc='pore.critical_temperature',
                       **params):
    r"""
    Uses Brock_Bird model to adjust surface tension from it's value at a given
    reference temperature to temperature of interest

    Parameters
    ----------
    To : float
        Reference temperature (K)

    sigma_o : float
        Surface tension at reference temperature (N/m)

    pore_T : string
        The dictionary key containing the temperature values (K)

    pore_Tc : string
        The dictionary key containing the critical temperature values (K)
    """
    Tc = phase[pore_Tc]
    Ti = phase[pore_T]
    Tro = To/Tc
    Tri = Ti/Tc
    value = sigma_o*(1-Tri)**(11/9)/(1-Tro)**(11/9)
    return value
