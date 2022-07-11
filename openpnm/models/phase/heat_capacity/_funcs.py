import numpy as np


__all__ = [
    'mixture_heat_capacity_XXX',
]


def mixture_heat_capacity_XXX(target):
    zs = [target['pore.mole_fraction.' + c.name] for c in target.components.values()]
    Cps = [c['pore.heat_capacity'] for c in target.components.values()]
    Cpmix = np.sum([zs[i]*Cps[i] for i in range(len(zs))], axis=0)
    return Cpmix
