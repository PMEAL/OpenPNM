import openpnm.models as mods

water = {
    'pore.molecular_weight': {
        'model': mods.misc.constant,
        'value': 0.01802,
        'regen_mode': 'deferred',
        },
    'pore.critical_pressure': {
        'model': mods.misc.constant,
        'value': 2.2064E7,
        'regen_mode': 'deferred',
        },
    'pore.critical_temperature': {
        'model': mods.misc.constant,
        'value': 647.1,
        'regen_mode': 'deferred',
        },
    'pore.critical_volume': {
        'model': mods.misc.constant,
        'value': 0.003106,
        'regen_mode': 'deferred',
        },
    'pore.contact_angle': {
        'model': mods.misc.constant,
        'value': 110.0,
        'regen_mode': 'deferred',
        },
    'pore.electrical_conductivity': {
        'model': mods.misc.constant,
        'value': 1e-15,
        'regen_mode': 'deferred',
        },
    'pore.diffusivity': {
        'model': mods.misc.constant,
        'value':  1e-9,
        'regen_mode': 'deferred',
        },
    'pore.density': {
        'model': mods.phases.density.water,
        'regen_mode': 'deferred',
        },
    'pore.molar_density': {
        'model': mods.phases.molar_density.standard,
        'regen_mode': 'deferred',
        },
    'pore.surface_tension': {
        'model': mods.phases.surface_tension.water,
        'regen_mode': 'deferred',
        },
    'pore.thermal_conductivity': {
        'model': mods.phases.thermal_conductivity.water,
        'regen_mode': 'deferred',
        },
    'pore.vapor_pressure': {
        'model': mods.phases.vapor_pressure.antoine,
        'A': 8.088,
        'B': 1750.71,
        'C': 236.191,
        'regen_mode': 'deferred',
        },
    'pore.viscosity': {
        'model': mods.phases.viscosity.water,
        'regen_mode': 'deferred',
        },
    }
