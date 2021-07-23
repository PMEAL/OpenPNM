import openpnm.models as mods


circles_and_rectangles = {
    'pore.seed': {
        'model': mods.misc.random,
        'element': 'pore',
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
    'pore.volume': {
        'model': mods.geometry.pore_volume.circle,
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
    'throat.length': {
        'model': mods.geometry.throat_length.circles_and_rectangles,
        'pore_diameter': 'pore.diameter',
        'throat_diameter': 'throat.diameter',
        'regen_mode': 'deferred',
        },
    'throat.volume': {
        'model': mods.geometry.throat_volume.rectangle,
        'throat_diameter': 'throat.diameter',
        'throat_length': 'throat.length',
        'regen_mode': 'deferred',
        },
    'throat.cross_sectional_area': {
        'model': mods.geometry.throat_cross_sectional_area.rectangle,
        'throat_diameter': 'throat.diameter',
        'regen_mode': 'deferred',
        },
    'throat.diffusive_size_factors': {
        'model': mods.geometry.diffusive_size_factors.circles_and_rectangles,
        'pore_diameter': 'pore.diameter',
        'throat_diameter': 'throat.diameter',
        'regen_mode': 'deferred',
        },
    'throat.hydraulic_size_factors': {
        'model': mods.geometry.hydraulic_size_factors.circles_and_rectangles,
        'pore_diameter': 'pore.diameter',
        'throat_diameter': 'throat.diameter',
        'regen_mode': 'deferred',
        },
    }

