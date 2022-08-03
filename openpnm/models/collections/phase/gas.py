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
        'model': mods.phase.heat_capacity.gas_pure_TRC,
    },
    'pore.heat_capacity': {
        'model': mods.phase.heat_capacity.gas_pure_TRC,
    },
    'pore.thermal_conductivity': {
        'model': mods.phase.thermal_conductivity.gas_pure_gismr,
    },
    'pore.viscosity': {
        'model': mods.phase.viscosity.gas_pure_gesmr,
    },
}


_gas_mixture = {
    'pore.density': {
        'model': mods.phase.density.ideal_gas,
    },
    'pore.heat_capacity': {
        'model': mods.phase.heat_capacity.gas_mixture_yweighted,
    },
    'pore.thermal_conductivity': {
        'model': mods.phase.thermal_conductivity.gas_mixture_whz,
    },
    'pore.viscosity': {
        'model': mods.phase.viscosity.gas_mixture_hz,
    },
    'pore.diffusivity': {
        'model': mods.phase.diffusivity.gas_mixture_ce,
    },
}
