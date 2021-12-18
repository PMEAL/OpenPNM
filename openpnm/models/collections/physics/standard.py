import openpnm.models.physics as mods

standard = {
    'throat.hydraulic_conductance': {
        'model': mods.hydraulic_conductance.generic_hydraulic,
        'regen_mode': 'deferred',
        },
    'throat.diffusive_conductance': {
        'model': mods.diffusive_conductance.generic_diffusive,
        'regen_mode': 'deferred',
        },
    'throat.entry_pressure': {
        'model': mods.capillary_pressure.washburn,
        'regen_mode': 'deferred',
        },
    }
