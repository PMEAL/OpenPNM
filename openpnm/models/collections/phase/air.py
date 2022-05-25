import openpnm.models as mods
from openpnm.utils import get_model_collection


def air(regen_mode='deferred', domain=None):
    return get_model_collection(collection=_air,
                                regen_mode=regen_mode,
                                domain=domain)


_air = {
    'pore.molecular_weight': {
        'model': mods.misc.constant,
        'value': 0.0291,
    },
    'pore.critical_pressure': {
        'model': mods.misc.constant,
        'value': 3.786E6,
    },
    'pore.critical_temperature': {
        'model': mods.misc.constant,
        'value': 132.5,
    },
    'pore.critical_volume': {
        'model': mods.misc.constant,
        'value': 0.002917,
    },
    'pore.contact_angle': {
        'model': mods.misc.constant,
        'value': 180.0,
    },
    'pore.surface_tension': {
        'model': mods.misc.constant,
        'value': 0.072,
    },
    'pore.molar_density': {
        'model': mods.phase.molar_density.ideal_gas,
    },
    'pore.diffusivity': {
        'model': mods.phase.diffusivity.fuller,
        'MA': 0.032,
        'MB': 0.028,
        'vA': 16.6,
        'vB': 17.9,
    },
    'pore.thermal_conductivity': {
        'model': mods.misc.polynomial,
        'prop': 'pore.temperature',
        'a': [0.00422791, 0.0000789606, -1.56383E-08],
    },
    'pore.electrical_conductivity': {
        'model': mods.misc.constant,
        'value': 1e-15,
    },
    'pore.viscosity': {
        'model': mods.misc.polynomial,
        'prop': 'pore.temperature',
        'a': [0.00000182082, 6.51815E-08, -3.48553E-11, 1.11409E-14],
    },
}
