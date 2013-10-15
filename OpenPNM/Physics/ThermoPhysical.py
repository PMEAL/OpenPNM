
"""
module ThermoPhysics
===============================================================================

"""

import OpenPNM
import scipy as sp

def Antoine(network,
            A=8.07131,
            B=1730.63,
            C=233.426):
    r"""
    ----
    """
    T = network.pore_properties['temperature']
    P = 10**(A - B/(C+T)) *(1.01325e5/760)
    return P
    
def HenrysLaw(network,Y):
    r"""
    ----
    """
    K = 4.31543175e9    
    X = Y/K
    return X