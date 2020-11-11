r"""
"""
import numpy as _np


def percolating_continua(target, phi_crit, tau,
                         volume_fraction='pore.volume_fraction',
                         bulk_property='pore.intrinsic_conductivity'):
    r'''
    Calculates the effective property of a continua using percolation theory

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    volume_fraction : string
        The dictionary key in the Phase object containing the volume fraction
        of the conducting component

    bulk_property : string
        The dictionary key in the Phase object containing the intrinsic
        property of the conducting component

    phi_crit : float
        The volume fraction below which percolation does NOT occur

    tau : float
        The exponent of the percolation relationship

    Returns
    -------
    sigma_eff : NumPy ndarray
        Array containing effective electrical conductivity values.

    Notes
    -----
    This model uses the following standard percolation relationship:

    .. math::

        \sigma_{effective}=\sigma_{bulk}(\phi - \phi_{critical})^\lambda

    '''
    sigma = target[bulk_property]
    phi = target[volume_fraction]
    diff_phi = _np.clip(phi - phi_crit, a_min=0, a_max=_np.inf)
    sigma_eff = sigma*(diff_phi)**tau
    return sigma_eff
