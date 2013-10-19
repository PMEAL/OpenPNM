
"""
module ContactAngle
===============================================================================

"""

import scipy as _sp

def set_as(theta=0,fluid2={}):
    r"""
    Set the contact of fluid relative to fluid2

    Parameters
    ----------
    theta : array_like
        The numerical value of contact angle in degrees

    fluid2 : OpenPNM Fluid Object
        The fluid against which this contact and applies
    """
    if fluid2:
        return {'contact_angle': {fluid2['name']: theta}}
    else:
        raise Exception('Defining contact angle requires specifying a second fluid')

