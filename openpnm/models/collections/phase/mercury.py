import openpnm.models as mods

mercury = {
    'pore.vapor_pressure': {
        'model': mods.phases.vapor_pressure.antoine,
        'A': 9.85767,
        'B': 3007.129,
        'C': -10.001,
        'regen_mode': 'deferred',
        },
    'pore.density': {
        'model': mods.misc.linear,
        'prop': 'pore.temperature',
        'b': 14280.9,
        'm': -2.47004,
        'regen_mode': 'deferred',
        },
    'pore.molar_density': {
        'model': mods.phases.molar_density.standard,
        'regen_mode': 'deferred',
        },
    'pore.surface_tension': {
        'model': mods.misc.linear,
        'prop': 'pore.temperature',
        'b': 0.56254,
        'm': -0.00028,
        'regen_mode': 'deferred',
        },
    'pore.thermal_conductivity': {
        'model': mods.misc.polynomial,
        'prop': 'pore.temperature',
        'a': [3.98691, 0.0170967, -0.0000063862],
        'regen_mode': 'deferred',
        },
    'pore.viscosity': {
        'model': mods.misc.polynomial,
        'prop': 'pore.temperature',
        'a': [0.00355837, -0.0000100131, 1.23684E-08, -5.1684E-12],
        'regen_mode': 'deferred',
        },
    }
