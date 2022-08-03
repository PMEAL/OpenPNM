import openpnm.models as mods
from openpnm.utils import get_model_collection


def water(regen_mode=None, domain=None):
    return get_model_collection(collection=_water,
                                regen_mode=regen_mode,
                                domain=domain)


_water = {
    'pore.contact_angle': {
        'model': mods.misc.constant,
        'value': 110.0,
    },
    'pore.density': {
        'model': mods.phase.density.water_correlation,
    },
    'pore.molar_density': {
        'model': mods.phase.density.mass_to_molar,
    },
    'pore.surface_tension': {
        'model': mods.phase.surface_tension.water_correlation,
    },
    'pore.thermal_conductivity': {
        'model': mods.phase.thermal_conductivity.water_correlation,
    },
    'pore.vapor_pressure': {
        'model': mods.phase.vapor_pressure.liquid_pure_antoine,
    },
    'pore.viscosity': {
        'model': mods.phase.viscosity.water_correlation,
    },
}
