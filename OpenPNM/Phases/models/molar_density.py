r"""
===============================================================================
Submodule -- molar_density
===============================================================================

"""
import scipy as sp

def standard(phase,**kwargs):
    r"""
    Calculates the molar density from the molecular with and mass density
    """
    MW = phase['pore.molecular_weight']
    rho = phase['pore.density']
    value = rho/MW
    return value

def ideal_gas(phase,**kwargs):
    r"""
    Uses ideal gas law to calculate the molar density of an ideal gas
 
    Parameters
    ----------
    P, T, MW: float, array_like
        P pressure of the gas in [Pa]
        T temperature of the gas in [K]
        MW molecular weight of the gas in [kg/kmole]
            
    Returns
    -------
    rho, the density in [mol/m3]
    
    """
   
    P = phase['pore.pressure']
    T = phase['pore.temperature']
    R = 8.31447
    value = P/(R*T)
    return value

def vanderwaals(phase,P,T,Pc,Tc,MW,**kwargs):
    r"""
    Uses Van der Waals equation of state to calculate the density of a real gas
 
    Parameters
    ----------
    P, T, Pc, Tc, MW: float, array_like
        P pressure of the gas in [Pa]
        T temperature of the gas in [K]
        Pc critical pressure of the gas in [Pa]
        T critical temperature of the gas in [K]
        MW molecular weight of the gas in [kg/mol]
            
    Returns
    -------
    rho, the density in [mol/m3]
    
    """
    
    P = phase['pore.pressure']
    T = phase['pore.temperature']
    Pc = phase['pore.criticalpressure']
    Tc = phase['pore.criticaltemperature']
    R = 8.314
    a = 27*(R**2)*(Tc**2)/(64*Pc); b = R*Tc/(8*Pc)
    a0 = 1; a1 = -1/b; a2 = (R*T+b*P)/(a*b); a3 = -P/(a*b)
    density = sp.roots([a0, a1, a2, a3])
    value = sp.real(density[2])
    return value