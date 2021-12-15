import openpnm.models as mods


stick_and_ball = {
    'pore.seed': {
        'model': mods.misc.random,
        'element':'pore',
        'num_range': [0.2, 0.7],
        'seed': None,
        'regen_mode': 'deferred',
        },
    'pore.max_size': {
        'model': mods.geometry.pore_size.largest_sphere,
        'iters': 10,
        'regen_mode': 'deferred',
        },
    'pore.diameter': {
        'model': mods.misc.product,
        'prop1': 'pore.max_size',
        'prop2': 'pore.seed',
        'regen_mode': 'deferred',
        },
    'pore.area': {
        'model': mods.geometry.pore_cross_sectional_area.sphere,
        'pore_diameter': 'pore.diameter',
        'regen_mode': 'deferred',
        },
    'pore.volume': {
        'model': mods.geometry.pore_volume.sphere,
        'pore_diameter': 'pore.diameter',
        'regen_mode': 'deferred',
        },
    'throat.max_size': {
        'model': mods.misc.from_neighbor_pores,
        'mode': 'min',
        'prop': 'pore.diameter',
        'regen_mode': 'deferred',
        },
    'throat.diameter': {
        'model': mods.misc.scaled,
        'factor': 0.5,
        'prop': 'throat.max_size',
        'regen_mode': 'deferred',
        },
    'throat.endpoints': {
        'model': mods.geometry.throat_endpoints.spherical_pores,
        'pore_diameter': 'pore.diameter',
        'throat_diameter': 'throat.diameter',
        'regen_mode': 'deferred',
        },
    'throat.length': {
        'model': mods.geometry.throat_length.piecewise,
        'throat_endpoints': 'throat.endpoints',
        'regen_mode': 'deferred',
        },
    'throat.surface_area': {
        'model': mods.geometry.throat_surface_area.cylinder,
        'throat_diameter': 'throat.diameter',
        'throat_length': 'throat.length',
        'regen_mode': 'deferred',
        },
    'throat.volume': {
        'model': mods.geometry.throat_volume.cylinder,
        'throat_diameter': 'throat.diameter',
        'throat_length': 'throat.length',
        'regen_mode': 'deferred',
        },
    'throat.area': {
        'model': mods.geometry.throat_cross_sectional_area.cylinder,
        'throat_diameter': 'throat.diameter',
        'regen_mode': 'deferred',
        },
    'throat.conduit_lengths': {
        'model': mods.geometry.throat_length.conduit_lengths,
        'throat_endpoints': 'throat.endpoints',
        'throat_length': 'throat.length',
        'regen_mode': 'deferred',
        },
     }
