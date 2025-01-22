import openpnm.pnmlib as pnm
import numpy as np
import matplotlib.pyplot as plt
from openpnm.pnmlib.inspect import tree


# Create an empty project with just a plain dict
prj = {}

# Generate network topology
Nx, Ny, Nz = 2, 2, 1
pn = pnm.generators.cubic(
    shape=[Nx, Ny, Nz],
    node_prefix='pore',
    edge_prefix='throat',
    )

pn['pore.left'] = pn['pore.coords'][:, 0] < 1
pn['pore.right'] = pn['pore.coords'][:, 0] > (Nx - 1)

pnm.core.add_network(prj, pn)


onetime_models = {
    'network/pore.seed': {
        'model': pnm.models.geometry.random_seeds,
        'num_range': [0.2, 0.8],
    },
    'network/pore.volume@left': {
        'model': pnm.models.geometry.constant,
        'value': 3.3,
    },
    'network/pore.volume@right': {
        'model': pnm.models.geometry.constant,
        'value': 1.1,
    },
}

geo_models1 = {
    'pore.size3': {
        'model': pnm.models.geometry.product,
        'props': ['pore.seed', 'pore.seed'],
    },
}

geo_models2 = {
    'network/pore.seed': {
        'model': pnm.models.geometry.random_seeds,
        'num_range': [0.2, 0.8],
    },
    'network/pore.size2@left': {
        'model': pnm.models.geometry.product,
        'props': ['network/pore.seed', 'network/pore.seed'],
    },
}


pnm.core.set_data(prj, 'param.T', 298.0)
pnm.core.set_data(prj, 'param.P', 101325)

pnm.core.set_data(prj, 'network/throat.size', 1.0)

pnm.models.apply_models(prj, onetime_models)
pnm.models.apply_models(prj, geo_models1, group='network')
pnm.models.apply_models(prj, geo_models2)


pnm.core.set_data(prj, 'phase1/pore.temperature', 1.0)
pnm.core.set_data(prj, 'phase1/comp1/pore.temperature', 1.0)

pnm.inspect.tree(prj)
pnm.inspect.info(pnm.core.get_data(prj, 'phase1/*'))

# phase.update(pnm.core.get_data(pn, '*.all'))

# phase_models = {
#     'pore.viscosity': {
#         'model': pnm.models.geometry.constant,
#         'value': 0.001,
#     },
#     'throat.viscosity': {
#         'model': pnm.models.geometry.constant,
#         'value': 0.001,
#     },
# }

# phase = pnm.models.apply_models(phase, phase_models)

# phys_models = {
#     'throat.hydraulic_conductance': {
#         'model': pnm.models.physics.hydraulic_conductance,
#         'network': pn,
#         'viscosity': 'viscosity',
#         'diameter': 'size',
#     },
# }

# phase = pnm.core.create_phase(prj)
# phase.update(pnm.core.get_data(pn, '*.all'))

# assert pn is prj['network']
# assert phase is prj['phase_02']

# pnm.core.set_data(prj['phase_02'], 'pore.viscosity', 0.001)
# pnm.core.set_data(prj, 'phase_02/throat.viscosity', 0.001)

# phys_models = {
#     'phase_02/throat.hydraulic_conductance': {
#         'model': pnm.models.physics.hydraulic_conductance2,
#         'pore_viscosity': 'phase_02/pore.viscosity',
#         'throat_viscosity': 'phase_02/throat.viscosity',
#         'pore_diameter': 'network/pore.size',
#         'throat_diameter': 'network/throat.size',
#     },
# }

# pnm.models.apply_models(prj, phys_models)

# A = pnm.simulations.build_A(prj, 'phase_02/throat.hydraulic_conductance')
# b = pnm.simulations.build_b(pn)
# Ps = pnm.core.get_pores(pn, 'left')
# A, b = pnm.simulations.set_value_bc(A, b, values=1.0, locs=Ps)
# Ps = pnm.core.get_pores(pn, 'right')
# A, b = pnm.simulations.set_value_bc(A, b, values=0.0, locs=Ps)
# x = pnm.simulations.solve(A, b)
# phase['pore.pressure'] = x
# # plt.imshow(x.reshape([Nx, Ny]))

# phase2 = pnm.core.create_phase(prj)
# phase2.update(pnm.core.get_data(pn, '*.all'))
# phase3 = pnm.core.create_phase(prj['phase_02'])
# phase3.update(pnm.core.get_data(pn, '*.all'))

# pnm.inspect.tree(prj)
# pnm.inspect.tree(prj, hide=['data'])

# pnm.inspect.info(pn)








