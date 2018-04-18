import openpnm as op
import scipy as sp

ws = op.core.Workspace()
ws.settings['local_data'] = True

# NETWORK
sp.random.seed(17)
nx, ny, nz = 30, 20, 1
pn = op.network.Cubic(shape=[nx, ny, nz], spacing=1e-4, name='pn11')

# GEOMETRIES
geom = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)

# PHASES
water = op.phases.Water(network=pn)

# PHYSICS
phys_water = op.physics.GenericPhysics(network=pn, phase=water, geometry=geom)

water['throat.viscosity'] = water['pore.viscosity'][0]
mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys_water.add_model(propname='throat.hydraulic_conductance',
                     model=mod, viscosity='throat.viscosity')

geom['pore.area'] = sp.pi*(geom['pore.diameter']**2)/4.0
mod2 = op.models.physics.diffusive_conductance.bulk_diffusion
phys_water.add_model(propname='throat.diffusive_conductance',
                     model=mod2, diffusivity='pore.diffusivity')

phys_water.regenerate_models()

inlet = pn.pores('back')  # pore inlet
outlet = pn.pores('front')  # pore outlet

inlet2 = pn.pores('left')  # pore inlet2
outlet2 = pn.pores('right')  # pore outlet2

# ALGORITHMS
alg1 = op.algorithms.StokesFlow(network=pn, phase=water)
alg1.set_BC(pores=inlet, bctype='dirichlet', bcvalues=10)
alg1.set_BC(pores=outlet, bctype='dirichlet', bcvalues=0)
alg1.run()
water['pore.pressure'] = alg1['pore.pressure']

alg2 = op.algorithms.TransientAdvectionDiffusion(network=pn, phase=water,
                                                 t_scheme='cranknicolson')
alg2.set_IC(0)
alg2.set_BC(pores=inlet, bctype='dirichlet', bcvalues=2)
alg2.set_BC(pores=outlet, bctype='dirichlet', bcvalues=0)
alg2.run()

alg3 = op.algorithms.AdvectionDiffusion(network=pn, phase=water)
alg3.set_BC(pores=inlet, bctype='dirichlet', bcvalues=2)
alg3.set_BC(pores=outlet, bctype='dirichlet', bcvalues=0)
alg3.run()

test = alg3['pore.mole_fraction'] - alg2['pore.mole_fraction']
