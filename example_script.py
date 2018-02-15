import openpnm as op
import scipy as sp
ws = op.core.Workspace()
ws.settings['local_data'] = True

sp.random.seed(0)
pn = op.network.Cubic(shape=[5, 5, 5], spacing=0.0001, name='pn11')
pn.add_boundary_pores()
sim = pn.simulation

geom = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts,
                                settings={'test': 1})

air = op.phases.Air(network=pn)
water = op.phases.Water(network=pn)
water['throat.viscosity'] = water['pore.viscosity'][0]

mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys_water = op.physics.GenericPhysics(network=pn, phase=water, geometry=geom)
phys_water.add_model(propname='throat.conductance',
                     model=mod, viscosity='throat.viscosity')
phys_water.regenerate_models()

alg = op.algorithms.FickianDiffusion(network=pn, phase=water)
alg.setup(conductance='throat.conductance', quantity='pore.mole_fraction')
alg.set_BC(pores=pn.pores('top'), bctype='dirichlet', bcvalues=0.5)
alg.set_BC(pores=pn.pores('bottom'), bctype='dirichlet', bcvalues=0.0)
alg['pore.mole_fraction'] = 0

rxn = op.algorithms.GenericReaction(network=pn, pores=[70, 71])
rxn['pore.k'] = 1e-1
rxn['pore.alpha'] = 1
rxn.add_model(propname='pore.rxn_rate',
              model=op.algorithms.models.standard_kinetics,
              quantity='pore.mole_fraction',
              prefactor='pore.k', exponent='pore.alpha',
              regen_mode='deferred')
rxn.settings['rate_model'] = 'pore.rxn_rate'
alg.set_source(source=rxn)

rxn2 = op.algorithms.GenericReaction(network=pn, pores=[50, 51])
rxn2.models = rxn.models.copy()
rxn2.settings['rate_model'] = 'pore.rxn_rate'
rxn2['pore.k'] = 1e-5
rxn2['pore.alpha'] = 2
alg.set_source(source=rxn2)

alg.run()

# op.io.VTK.save(simulation=pn.simulation, phases=[water])
