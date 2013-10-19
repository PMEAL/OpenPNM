
"""
module SurfaceTension
===============================================================================

"""

import scipy as _sp

def set_as(sigma=0.072,fluid2={}):
    r"""
    Set the surface tension of the fluid relative to fluid2

    Parameters
    ----------
    sigma : array_like
        The numerical value of surface tension

    fluid2 : OpenPNM Fluid Object
        The fluid against which this surface tension applies
    """
    if fluid2:
        return {'surface_tension': {fluid2['name']: sigma}}
    else:
        raise Exception('Defining surface tension requires supplying a second fluid')

