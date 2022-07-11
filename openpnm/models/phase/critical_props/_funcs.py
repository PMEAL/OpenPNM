import numpy as np


__all__ = [
    "liquid_mixture_Vc_XXX",
    "liquid_mixture_Tc_XXX",
    "liquid_mixture_acentric_factor_XXX",
]


def liquid_mixture_Vc_XXX(target):
    xs = [target['pore.mole_fraction.' + c.name] for c in target.components.values()]
    Vcs = [c['param.critical_volume'] for c in target.components.values()]
    N = len(xs)  # Number of components
    Vm1 = np.sum([xs[i]*Vcs[i] for i in range(N)], axis=0)
    Vm2 = np.sum([xs[i]*(Vcs[i])**(2/3) for i in range(N)], axis=0)
    Vm3 = np.sum([xs[i]*(Vcs[i])**(1/3) for i in range(N)], axis=0)
    Vm = 0.25*(Vm1 + 3*(Vm2*Vm3))
    return Vm


def liquid_mixture_Tc_XXX(target, Vc='pore.critical_volume'):
    xs = [target['pore.mole_fraction.' + c.name] for c in target.components.values()]
    Tcs = [c['param.critical_temperature'] for c in target.components.values()]
    Vcs = [c['param.critical_volume'] for c in target.components.values()]
    Vm = target[Vc]
    N = len(xs)  # Number of components
    num = np.zeros_like(xs[0])
    for i in range(N):
        for j in range(N):
            VT = (Vcs[i]*Tcs[i]*Vcs[j]*Tcs[j])**0.5
            num += xs[i]*xs[j]*VT
    Tcm = num/Vm
    return Tcm


def liquid_mixture_acentric_factor_XXX(target):
    xs = [target['pore.mole_fraction.' + c.name] for c in target.components.values()]
    omegas = [c['param.acentric_factor'] for c in target.components.values()]
    omega = np.sum([omegas[i]*xs[i] for i in range(len(xs))], axis=0)
    return omega
