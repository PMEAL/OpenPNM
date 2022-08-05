import openpnm.models.physics as mods


standard = {
    'throat.hydraulic_conductance': {
        'model': mods.hydraulic_conductance.generic_hydraulic,
    },
    'throat.diffusive_conductance': {
        'model': mods.diffusive_conductance.generic_diffusive,
    },
    'throat.entry_pressure': {
        'model': mods.capillary_pressure.washburn,
    },
    'throat.thermal_conductance': {
        'model': mods.thermal_conductance.generic_thermal,
    },
    'throat.electrical_conductance': {
        'model': mods.electrical_conductance.generic_electrical,
    },
    'throat.ad_dif_conductance': {
        'model': mods.ad_dif_conductance.ad_dif,
    },
}
