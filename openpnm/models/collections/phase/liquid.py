import openpnm.models as mods
from openpnm.utils import get_model_collection


def standard_liquid(regen_mode=None, domain=None):
    return get_model_collection(collection=_pure_liquid,
                                regen_mode=regen_mode,
                                domain=domain)


def standard_liquid_mixture(regen_mode=None, domain=None):
    return get_model_collection(collection=_liquid_mixture,
                                regen_mode=regen_mode,
                                domain=domain)


_pure_liquid = {
    'pore.density': {
        'model': mods.phase.density.liquid_pure_COSTALD,
    },
    'pore.heat_capacity_gas': {
        'model': mods.phase.heat_capacity.gas_pure_TRC,
    },
    'pore.heat_capacity': {
        'model': mods.phase.heat_capacity.liquid_pure_rp,
    },
    'pore.thermal_conductivity': {
        'model': mods.phase.thermal_conductivity.liquid_pure_gismr,
    },
    'pore.viscosity': {
        'model': mods.phase.viscosity.liquid_pure_ls,
    },
    'pore.vapor_pressure': {
        'model': mods.phase.vapor_pressure.liquid_pure_antoine,
    },
}


_liquid_mixture = {
    'pore.density': {
        'model': mods.phase.density.liquid_mixture_COSTALD,
    },
    'pore.heat_capacity': {
        'model': mods.phase.heat_capacity.liquid_mixture_xweighted,
    },
    'pore.thermal_conductivity': {
        'model': mods.phase.thermal_conductivity.liquid_mixture_DIPPR9H,
    },
    'pore.viscosity': {
        'model': mods.phase.viscosity.liquid_mixture_xweighted,
    },
}
