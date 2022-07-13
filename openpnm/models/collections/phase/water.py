import openpnm.models as mods
from openpnm.utils import get_model_collection


def water(regen_mode=None, domain=None):
    return get_model_collection(collection=_water,
                                regen_mode=regen_mode,
                                domain=domain)


_water = {
    'pore.molecular_weight': {
        'model': mods.misc.constant,
        'value': 0.01802,
    },
    'pore.critical_pressure': {
        'model': mods.misc.constant,
        'value': 2.2064E7,
    },
    'pore.critical_temperature': {
        'model': mods.misc.constant,
        'value': 647.1,
    },
    'pore.critical_volume': {
        'model': mods.misc.constant,
        'value': 0.003106,
    },
    'pore.contact_angle': {
        'model': mods.misc.constant,
        'value': 110.0,
    },
    'pore.electrical_conductivity': {
        'model': mods.misc.constant,
        'value': 1e-15,
    },
    'pore.diffusivity': {
        'model': mods.misc.constant,
        'value': 1e-9,
    },
    'pore.density': {
        'model': mods.phase.density.water_correlation,
    },
    'pore.molar_density': {
        'model': mods.phase.molar_density.standard,
    },
    'pore.surface_tension': {
        'model': mods.phase.surface_tension.water_correlation,
    },
    'pore.thermal_conductivity': {
        'model': mods.phase.thermal_conductivity.water_correlation,
    },
    'pore.vapor_pressure': {
        'model': mods.phase.vapor_pressure.antoine,
        'A': 8.088,
        'B': 1750.71,
        'C': 236.191,
    },
    'pore.viscosity': {
        'model': mods.phase.viscosity.water_correlation,
    },
}
