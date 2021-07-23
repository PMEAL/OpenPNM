import openpnm.models as models


spheres_and_cylinders = {
	'pore.seed': {
		'model': mods.misc.random,
		'element': 'pore',
		'num_range': [0.2, 0.7],
		'seed': None,
		},
	'pore.max_size': {
		'model': mods.geometry.pore_size.largest_sphere,
		'iters': 10,
		},
	'pore.diameter': {
		'model': mods.misc.product,
		'prop1': 'pore.max_size',
        'prop2': 'pore.seed',
        },
	'pore.volume': {
		'model': mods.geometry.pore_volume.sphere,
		'pore_diameter': 'pore.diameter',
		}
	'throat.max_size': {
		'model': mods.misc.from_neighbor_pores,
		'mode': 'min',
        'prop': 'pore.diameter',
        },
	'throat.diameter': {
		'model': mods.misc.scaled,
		'factor': 0.5,
        'prop': 'throat.max_size',
        },
	'throat.length': {
		'model': mods.geometry.throat_length.spheres_and_cylinders,
		'pore_diameter': 'pore.diameter',
        'throat_diameter': 'throat.diameter',
        },
	'throat.cross_sectional_area': {
		'model': mods.geometry.throat_cross_sectional_area.cylinder,
		'throat_diameter': 'throat.diameter',
		},
	'throat.volume': {
		'model': mods.geometry.throat_volume.cylinder,
		'throat_diameter': 'throat.diameter',
        'throat_length': 'throat.length',
        },
	'throat.diffusive_size_factors': {
		'model': gmods.diffusive_size_factors.spheres_and_cylinders,
		'pore_diameter': 'pore.diameter',
        'throat_diameter': 'throat.diameter',
        },
	'throat.hydraulic_size_factors': {
		'model': gmods.hydraulic_size_factors.spheres_and_cylinders,
		'pore_diameter': 'pore.diameter',
		'throat_diameter': 'throat.diameter',
		},
	}