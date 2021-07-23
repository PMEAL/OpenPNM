import openpnm.models as mods

standard = {
    'throat.flow_shape_factors': {
        'model': mods.physics.flow_shape_factors.conical_frustum_and_stick,
        'regen_mode': 'deferred',
        },
    'throat.hydraulic_conductance': {
        'model': mods.physics.hydraulic_conductance.hagen_poiseuille,
        'regen_mode': 'deferred',
        },
    'throat.poisson_shape_factors': {
        'model': mods.physics.poisson_shape_factors.conical_frustum_and_stick,
        'regen_mode': 'deferred',
        },
    'throat.diffusive_conductance': {
        'model': mods.physics.diffusive_conductance.mixed_diffusion,
        'regen_mode': 'deferred',
        },
    'throat.ad_dif_conductance': {
        'model': mods.physics.ad_dif_conductance.ad_dif,
        'regen_mode': 'deferred',
        },
    'throat.entry_pressure': {
        'model': mods.physics.capillary_pressure.washburn,
        'regen_mode': 'deferred',
        },
    'throat.thermal_conductance': {
        'model': mods.physics.thermal_conductance.series_resistors,
        'regen_mode': 'deferred',
        },
    'throat.electrical_conductance': {
        'model': mods.physics.electrical_conductance.series_resistors,
        'regen_mode': 'deferred',
        },
    }