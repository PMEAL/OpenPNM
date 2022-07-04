import openpnm.models as mods
from openpnm.utils import get_model_collection


def gas_mixture(regen_mode=None, domain=None):
    return get_model_collection(collection=_gas_mixture,
                                regen_mode=regen_mode,
                                domain=domain)


_gas_mixture = {
    'pore.molecular_weight': {
        'model': mods.mixtures.mixture_molecular_weight,
    },
    'pore.gas_viscosity': {
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
