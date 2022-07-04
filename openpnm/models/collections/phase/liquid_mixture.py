import openpnm.models as mods
from openpnm.utils import get_model_collection


def liquid_mixture(regen_mode=None, domain=None):
    return get_model_collection(collection=_liquid_mixture,
                                regen_mode=regen_mode,
                                domain=domain)


_liquid_mixture = {
    'pore.molecular_weight': {
        'model': mods.mixtures.mixture_molecular_weight,
    },
    'pore.viscosity': {
        'model': mods.viscosity.liquid_mixture_viscosity,
    },
    'pore.critical_volume': {
        'model': mods.critical_props.liquid_mixture_critical_volume,
    },
    'pore.critical_temperature': {
        'model': mods.critical_props.liquid_mixture_critical_temperature,
    },
    'pore.acentric_factor': {
        'model': mods.critical_props.mixture_acentric_factor,
    },
    'pore.density': {
        'model': mods.density.liquid_mixture_density,
    },
    'pore.thermal_conductivity': {
        'model': mods.thermal_conductivity.liquid_mixture_thermal_conductivity,
    },
    'pore.heat_capacity': {
        'model': mods.heat_capacity.mixture_heat_capacity,
    },
}
