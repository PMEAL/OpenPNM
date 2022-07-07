import openpnm.models.phase as mods
from openpnm.utils import get_model_collection


def binary_gas_mixture(regen_mode=None, domain=None):
    return get_model_collection(collection=_binary_gas_mixture,
                                regen_mode=regen_mode,
                                domain=domain)


def gas_mixture(regen_mode=None, domain=None):
    return get_model_collection(collection=_gas_mixture,
                                regen_mode=regen_mode,
                                domain=domain)


def generic_gas(regen_mode=None, domain=None):
    return get_model_collection(collection=_generic_gas,
                                regen_mode=regen_mode,
                                domain=domain)


_generic_gas = {
    'pore.heat_capacity': {
        'model': mods.heat_capacity.gas_heat_capacity,
    },
    'pore.thermal_conductivity': {
        'model': mods.thermal_conductivity.gas_thermal_conductivity,
    },
    'pore.viscosity': {
        'model': mods.viscosity.gas_viscosity,
    },
}


_binary_gas_mixture = {
    'pore.molecular_weight': {
        'model': mods.mixtures.mixture_molecular_weight,
    },
    'pore.viscosity': {
        'model': mods.viscosity.gas_mixture_viscosity,
    },
    'pore.thermal_conductivity': {
        'model': mods.thermal_conductivity.gas_mixture_thermal_conductivity,
    },
    'pore.heat_capacity': {
        'model': mods.heat_capacity.mixture_heat_capacity,
    },
    'pore.LJ_epsilon': {
        'model': mods.diffusivity.gas_mixture_LJ_epsilon,
    },
    'pore.LJ_sigma': {
        'model': mods.diffusivity.gas_mixture_LJ_sigma,
    },
    'pore.LJ_omega': {
        'model': mods.diffusivity.gas_mixture_LJ_collision_integral,
    },
    'pore.diffusivity': {
        'model': mods.diffusivity.gas_mixture_diffusivity,
    },
}


_gas_mixture = {
    'pore.molecular_weight': {
        'model': mods.mixtures.mixture_molecular_weight,
    },
    'pore.viscosity': {
        'model': mods.viscosity.gas_mixture_viscosity,
    },
    'pore.thermal_conductivity': {
        'model': mods.thermal_conductivity.gas_mixture_thermal_conductivity,
    },
    'pore.heat_capacity': {
        'model': mods.heat_capacity.mixture_heat_capacity,
    },
}
