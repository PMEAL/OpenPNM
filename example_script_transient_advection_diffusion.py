import openpnm as op
import scipy as sp
import matplotlib.pyplot as plt

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
alg1.set_dirichlet_BC(pores=inlet, values=10)
alg1.set_dirichlet_BC(pores=outlet, values=0)
alg1.run()
water['pore.pressure'] = alg1['pore.pressure']

alg2 = op.algorithms.TransientAdvectionDiffusion(network=pn, phase=water,
                                                 t_scheme='cranknicolson')
alg2.set_IC(0)
alg2.set_dirichlet_BC(pores=inlet, values=2)
alg2.set_dirichlet_BC(pores=outlet, values=0)
alg2.run()

# PLOT
Z = sp.array([sp.reshape(alg2['pore.mole_fraction_initial'], (nx, ny)),
              sp.reshape(alg2['pore.mole_fraction0'], (nx, ny)),
              sp.reshape(alg2['pore.mole_fraction5'], (nx, ny)),
              sp.reshape(alg2['pore.mole_fraction10'], (nx, ny)),
              sp.reshape(alg2['pore.mole_fraction20'], (nx, ny)),
              sp.reshape(alg2['pore.mole_fraction40'], (nx, ny)),
              sp.reshape(alg2['pore.mole_fraction60'], (nx, ny)),
              sp.reshape(alg2['pore.mole_fraction_steady'], (nx, ny))])

fig, axes = plt.subplots(nrows=2, ncols=4)
i = 0
for ax in axes.flat:
    im = ax.imshow(Z[i].T, cmap='rainbow')
    ax.set_title(str(i))
    ax.set_ylabel('y')
    ax.set_xlabel('x')
    i += 1

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
fig.colorbar(im, cax=cbar_ax)

plt.suptitle('OpenPNM transient dispersion', fontsize=16)
plt.show()
