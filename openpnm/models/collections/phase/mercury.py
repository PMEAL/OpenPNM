import openpnm.models as mods
from openpnm.utils import get_model_collection


def mercury(regen_mode=None, domain=None):
    return get_model_collection(collection=_mercury,
                                regen_mode=regen_mode,
                                domain=domain)


_mercury = {
    'pore.vapor_pressure': {
        'model': mods.phase.vapor_pressure.antoine,
        'A': 9.85767,
        'B': 3007.129,
        'C': -10.001,
    },
    'throat.contact_angle': {
        'model': mods.misc.constant,
        'value': 140,
    },
    'pore.molecular_weight': {
        'model': mods.misc.constant,
        'value': 200.59,
    },
    'pore.density': {
        'model': mods.misc.linear,
        'prop': 'pore.temperature',
        'b': 14280.9,
        'm': -2.47004,
    },
    'pore.molar_density': {
        'model': mods.phase.molar_density.mass_to_molar,
    },
    'pore.surface_tension': {
        'model': mods.misc.linear,
        'prop': 'pore.temperature',
        'b': 0.56254,
        'm': -0.00028,
    },
    'pore.thermal_conductivity': {
        'model': mods.misc.polynomial,
        'prop': 'pore.temperature',
        'a': [3.98691, 0.0170967, -0.0000063862],
    },
    'pore.viscosity': {
        'model': mods.misc.polynomial,
        'prop': 'pore.temperature',
        'a': [0.00355837, -0.0000100131, 1.23684E-08, -5.1684E-12],
    },
}
