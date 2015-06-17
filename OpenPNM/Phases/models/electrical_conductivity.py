r"""
===============================================================================
Submodule -- electrical_conductivity
===============================================================================

"""
import scipy as _sp


def percolating_continua(phase,
                         phi_crit,
                         tau,
                         volume_fraction='pore.volume_fraction',
                         bulk_property='pore.intrinsic_conductivity',
                         **kwargs):
    r'''
    Calculates the effective property of a continua using percolation theory

    Parameters
    ----------
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

    Notes
    -----
    This model uses the following standard percolation relationship:

    .. math::

        \sigma_{effective}=\sigma_{bulk}(\phi - \phi_{critical})^\lambda

    '''
    sigma = phase[bulk_property]
    phi = phase[volume_fraction]
    diff_phi = _sp.clip(phi - phi_crit, a_min=0, a_max=_sp.inf)
    sigma_eff = sigma*(diff_phi)**tau
    return sigma_eff
