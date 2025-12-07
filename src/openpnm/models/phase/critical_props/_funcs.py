import numpy as np
from openpnm.models.phase import _phasedocs


__all__ = [

]


@_phasedocs
def liquid_mixture_Vc_XXX(
    phase,
    Vcs='pore.critical_volume.*',
):  # pragma: no cover
    r"""
    Calculates the critical volume of a liquid mixture using the correlation
    in [1]_

    Parameters
    ----------
    %(phase)s
    %(Vcs)s

    Returns
    -------

    References
    ----------
    .. [1] tbd

    """
    raise NotImplementedError("This function is not ready yet")
    xs = phase['pore.mole_fraction']
    Vcs = phase.get_comp_vals(Vcs)
    N = len(xs)  # Number of components
    Vm1 = np.sum([xs[i]*Vcs[i] for i in range(N)], axis=0)
    Vm2 = np.sum([xs[i]*(Vcs[i])**(2/3) for i in range(N)], axis=0)
    Vm3 = np.sum([xs[i]*(Vcs[i])**(1/3) for i in range(N)], axis=0)
    Vm = 0.25*(Vm1 + 3*(Vm2*Vm3))
    return Vm


@_phasedocs
def liquid_mixture_Tc_XXX(
    phase,
    Vm='pore.molar_volume',
    Vcs='pore.critical_volume.*',
    Tcs='pore.critical_temperature.*',
):  # pragma: no cover
    r"""
    Calculates the critical temperature of a liquid mixture using the
    correlation in [1]_

    Parameters
    ----------
    %(phase)s
    %(Vm)s
    %(Vcs)s
    %(Tcs)s

    Returns
    -------

    References
    ----------
    .. [1] tbd

    """
    raise NotImplementedError("This function is not ready yet")
    xs = phase['pore.mole_fraction']
    Tcs = phase.get_comp_vals(Tcs)
    Vcs = phase.get_comp_vals(Vcs)
    Vm = phase[Vm]
    N = len(xs)  # Number of components
    num = np.zeros_like(xs[0])
    for i in range(N):
        for j in range(N):
            VT = (Vcs[i]*Tcs[i]*Vcs[j]*Tcs[j])**0.5
            num += xs[i]*xs[j]*VT
    Tcm = num/Vm
    return Tcm


def liquid_mixture_acentric_factor_XXX(
    phase,
    omegas='param.acentric_factor.*',
):  # pragma: no cover
    r"""
    Calculates the accentric factor of a liquid mixture using the correlation
    in [1]_

    Parameters
    ----------
    %(phase)s
    %(omegas)s

    Returns
    -------

    References
    ----------
    .. [1] tbd

    """
    raise NotImplementedError("This function is not ready yet")
    xs = phase['pore.mole_fraction']
    omegas = phase.get_comp_vals(omegas)
    omega = np.sum([omegas[i]*xs[i] for i in range(len(xs))], axis=0)
    return omega
