r"""
===============================================================================
Submodule -- molar_density
===============================================================================

"""
import scipy as sp

def standard(phase,**kwargs):
    r"""
    Calculates the molar density of a fluid from the density and 
    molecular weight values

    """
    rho = phase['pore.density']
    MW = phase['pore.molecular_weight']
    value = rho/MW
    return value
