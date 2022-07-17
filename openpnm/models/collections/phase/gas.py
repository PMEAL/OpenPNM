import openpnm.models as mods
from openpnm.utils import get_model_collection


def standard_gas(regen_mode=None, domain=None):
    return get_model_collection(collection=_pure_gas,
                                regen_mode=regen_mode,
                                domain=domain)


def standard_gas_mixture(regen_mode=None, domain=None):
    return get_model_collection(collection=_gas_mixture,
                                regen_mode=regen_mode,
                                domain=domain)


_pure_gas = {
    'pore.density': {
        'model': mods.phase.density.ideal_gas,
    },
    'pore.heat_capacity_gas': {
        'model': mods.phase.heat_capacity.gas_pure,
    },
    'pore.heat_capacity': {
        'model': mods.phase.heat_capacity.gas_pure,
    },
    'pore.thermal_conductivity': {
        'model': mods.phase.thermal_conductivity.gas_pure,
    },
    'pore.viscosity': {
        'model': mods.phase.viscosity.gas_pure,
    },
}


_gas_mixture = {
    'pore.density': {
        'model': mods.phase.density.ideal_gas,
    },
    'pore.heat_capacity': {
        'model': mods.phase.heat_capacity.gas_mixture,
    },
    'pore.thermal_conductivity': {
        'model': mods.phase.thermal_conductivity.gas_mixture,
    },
    'pore.viscosity': {
        'model': mods.phase.viscosity.gas_mixture,
    },
    'pore.diffusivity': {
        'model': mods.phase.diffusivity.gas_mixture_chapman_enskog,
    },
}
