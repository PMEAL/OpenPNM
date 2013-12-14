
"""
module vapor_pressure
===============================================================================

"""
import scipy as sp

def constant(fluid,value=3000,**params):
    fluid.pore_conditions['vapor_pressure'] = value

def na(fluid,**params):
    fluid.pore_conditions['vapor_pressure'] = -1

def Antoine(fluid,A=8.07131,B=1730.63,C=233.426,**params):
    r"""
    Uses Antoine equation to estimate vapor pressure of a pure component

    Parameters
    ----------
    A, B, C :  float, array_like
            Antoine vapor pressure constants for pure compounds

    """
    T = fluid.pore_conditions['temperature']
    Pv = (10**(A-B/(C+T-273.15)))*1.333e2
    fluid.pore_conditions['vapor_pressure'] = Pv