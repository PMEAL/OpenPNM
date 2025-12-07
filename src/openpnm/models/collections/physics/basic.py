import openpnm.models.physics as mods


basic = {
    'throat.hydraulic_conductance': {
        'model': mods.hydraulic_conductance.generic_hydraulic,
    },
    'throat.diffusive_conductance': {
        'model': mods.diffusive_conductance.generic_diffusive,
    },
    'throat.entry_pressure': {
        'model': mods.capillary_pressure.washburn,
    },
}
