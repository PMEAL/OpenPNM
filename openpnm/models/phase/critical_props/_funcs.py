import numpy as np


__all__ = [

]


def liquid_mixture_Vc_XXX(
    target,
    Vcs='pore.critical_volume.*',
):
    raise NotImplementedError("This function is not ready yet")
    xs = target['pore.mole_fraction']
    Vcs = target.get_comp_vals(Vcs)
    N = len(xs)  # Number of components
    Vm1 = np.sum([xs[i]*Vcs[i] for i in range(N)], axis=0)
    Vm2 = np.sum([xs[i]*(Vcs[i])**(2/3) for i in range(N)], axis=0)
    Vm3 = np.sum([xs[i]*(Vcs[i])**(1/3) for i in range(N)], axis=0)
    Vm = 0.25*(Vm1 + 3*(Vm2*Vm3))
    return Vm


def liquid_mixture_Tc_XXX(
    target,
    Vm='pore.molar_volume',
    Vcs='pore.critical_volume.*',
    Tcs='pore.critical_temperature.*',
):
    raise NotImplementedError("This function is not ready yet")
    xs = target['pore.mole_fraction']
    Tcs = target.get_comp_vals(Tcs)
    Vcs = target.get_comp_vals(Vcs)
    Vm = target[Vm]
    N = len(xs)  # Number of components
    num = np.zeros_like(xs[0])
    for i in range(N):
        for j in range(N):
            VT = (Vcs[i]*Tcs[i]*Vcs[j]*Tcs[j])**0.5
            num += xs[i]*xs[j]*VT
    Tcm = num/Vm
    return Tcm


def liquid_mixture_acentric_factor_XXX(
    target,
    omegas='param.acentric_factor.*',
):
    raise NotImplementedError("This function is not ready yet")
    xs = target['pore.mole_fraction']
    omegas = target.get_comp_vals(omegas)
    omega = np.sum([omegas[i]*xs[i] for i in range(len(xs))], axis=0)
    return omega
