import openpnm.models as mods

standard = {
    'throat.hydraulic_conductance': {
        'model': mods.physics.hydraulic_conductance.generic_hydraulic,
        'regen_mode': 'deferred',
        },
    'throat.diffusive_conductance': {
        'model': mods.physics.diffusive_conductance.generic_diffusive,
        'regen_mode': 'deferred',
        },
    'throat.entry_pressure': {
        'model': mods.physics.capillary_pressure.washburn,
        'regen_mode': 'deferred',
        },
    }
