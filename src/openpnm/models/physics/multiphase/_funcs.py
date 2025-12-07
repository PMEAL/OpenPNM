import numpy as np
from openpnm.models import _doctxt


__all__ = ["conduit_conductance", "late_filling"]


@_doctxt
def conduit_conductance(phase, throat_conductance,
                        throat_occupancy='throat.occupancy',
                        pore_occupancy='pore.occupancy',
                        mode='strict', factor=1e-6):
    r"""
    Determines the conductance of a pore-throat-pore conduit based on the
    invaded state of each element.

    Parameters
    ----------
    %(phase)s
    throat_conductance : str
        %(dict_blurb)s throat conductance for transport
    pore_occupancy : str
        %(dict_blurb)s pore occupancy of the phase associated with ``phase``.
        An occupancy of 1 means the pore is completely filled with the phase
        and it fully conducts.
    throat_occupancy : str
        %(dict_blurb)s throat occupancy of the phase associated with
        ``phase``. An occupancy of 1 means the pore is completely filled
        with the phase and it fully conducts.
    mode : str
        How agressive the method should be when determining if a conduit is
        closed. Options are:

        ========= ============================================================
        Mode      Description
        ========= ============================================================
        strict    If any pore or throat in the conduit is unoccupied by
                  the given phase, the conduit is closed.
        medium    If either the throat or both pores are unoccupied, the
                  conduit is closed
        loose     Only close the conduit if the throat is unoccupied
        ========= ============================================================

    factor : float (default is 1e-6)
        The factor which becomes multiplied to the original conduit's
        conductance to severely limit transport, but not set it to zero.

    Returns
    -------
    %(return_arr)s adjusted conductance values

    """
    network = phase.network
    Tinv = phase[throat_occupancy] < 0.5
    P12 = network['throat.conns']
    Pinv = phase[pore_occupancy][P12] < 0.5
    if mode == 'loose':
        mask = Tinv
    elif mode == 'medium':
        mask = Tinv + np.all(Pinv, axis=1)
    elif mode == 'strict':
        mask = Tinv + np.any(Pinv, axis=1)
    else:
        raise Exception('Unrecongnized mode '+mode)
    value = phase[throat_conductance].copy()
    value[mask] = value[mask]*factor
    return value


@_doctxt
def late_filling(phase, pressure='pore.pressure',
                 Pc_star='pore.pc_star',
                 Swp_star=0.2, eta=3):
    r"""
    Calculates the fraction of a pore or throat filled with invading fluid
    based on the capillary pressure in the invading phase. The invading phase
    volume is calculated from:

        .. math::
            S_{nwp} = 1 - S_{wp}^{*} (P^{*}/P_{c})^{\eta}

    Parameters
    ----------
    pressure : str
        %(dict_blurb)s capillary pressure in the non-wetting phase (Pc > 0).
    Pc_star : str
        %(dict_blurb)s minimum pressure required to create an interface
        within the pore body or throat.  Typically this would be calculated
        using the Washburn equation.
    Swp_star : float
        The residual wetting phase in an invaded pore or throat at a pressure
        of ``pc_star``.
    eta : float
        Exponent controlling the rate at which wetting phase is displaced with
        increasing pressure.

    Returns
    -------
    %(return_arr)s fraction of each pore or throat that would be filled with
    non-wetting phase at the given phase pressure. This does not account
    for whether or not the element is actually invaded, which requires a
    percolation algorithm of some sort.

    """
    pc_star = phase[Pc_star]
    Pc = phase[pressure]
    # Remove any 0's from the Pc array to prevent numpy div by 0 warning
    Pc = np.maximum(Pc, 1e-9)
    Swp = Swp_star*((pc_star/Pc)**eta)
    values = np.clip(1 - Swp, 0.0, 1.0)
    return values
