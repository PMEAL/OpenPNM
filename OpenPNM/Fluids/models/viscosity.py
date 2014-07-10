r"""
===============================================================================
Submodule -- viscosity
===============================================================================

Models for predicting fluid viscosity

"""
import scipy as _sp

def reynolds(fluid,uo,b,**kwargs):
    r"""
    Uses Reynolds method for the temperature dependance of shear viscosity

    Parameters
    ----------
    u0, b :  float, array_like
            Coefficients of Reynolds method

    """
    T = fluid['pore.temperature']
    value = uo*_sp.exp(-1*b*T)
    return value

def chung(fluid,Vc,MW,acentric,kappa,dipole,**kwargs):
    r"""
    Uses Chung et al. model to estimate viscosity for gases with low pressure 
    (not near the critical pressure) from first principles at conditions of 
    interest

    Parameters
    ----------
    Vc :  float, array_like
        Critical volume of the component (m3/mol)
    MW : float, array_like
        Molecular weight of the component (kg/mol)
    acentric : float, array_like
        Acentric factor of the component
    kappa : float, array_like
        Special correction for highly polar substances
    dipole :float, array_like
        Dipole moment (C.m)

    """
    T = fluid['pore.temperature']
    Tc = fluid['fluid.Tc']
    Tr= T/Tc
    Tstar = 1.2593*Tr
    A = 1.161415
    B = 0.14874
    C = 0.52487
    D = 0.77320
    E = 2.16178
    F = 2.43787
    sigma = (A*(Tstar)**(-B)) + C*(_sp.exp(-D*Tstar)) + E*(_sp.exp(-F*Tstar))
    dipole_r = 131.3*(dipole*2.997e29)/((Vc*1e6)*Tc)**0.5
    f = 1-0.2756*acentric + 0.059035*(dipole_r**4) + kappa
    value = 40.785*f*(MW*1e3*T)**0.5/((Vc*1e6)**(2/3)*sigma)*1e-7
    return value

