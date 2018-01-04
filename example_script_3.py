import openpnm as op
import scipy as sp
ws = op.core.Workspace()
ws.settings['local_data'] = True

sp.random.seed(0)
pn = op.network.Cubic(shape=[15, 15, 15], spacing=0.0001, name='pn')
# pn.add_boundaries()

geom = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)

water = op.phases.Water(network=pn)
water['throat.viscosity'] = water['pore.viscosity'][0]
water.add_model(propname='throat.conductance',
                model=op.physics.models.hydraulic_conductance.hagen_poiseuille,
                viscosity='throat.viscosity')
water.regenerate_models()

alg = op.algorithms.FickianDiffusion(network=pn, phase=water)
alg.setup(conductance='throat.conductance', quantity='pore.mole_fraction')
alg.set_boundary_conditions(pores=pn.pores('top'), bctype='dirichlet',
                            bcvalues=0.5)
# alg.set_boundary_conditions(pores=pn.pores('bottom'), bctype='dirichlet',
#                             bcvalues=20)

alg.build_A()
alg.build_b()

Ps = [888, 889]
rxn = op.algorithms.GenericReaction(network=pn, algorithm=alg)
rxn['pore.A'] = 0.0
rxn['pore.A'][Ps] = 1e-5
rxn['pore.b'] = 2
rxn.add_model(propname='pore.rate',
              model=op.algorithms.models.standard_kinetics,
              quantity='pore.mole_fraction',
              prefactor='pore.A', exponent='pore.b')
rxn.setup(quantity='pore.mole_fraction')
rxn.solve()

# x = alg['pore.pressure']
# water['pore.pressure'] = x

# op.io.VTK.save(simulation=pn.simulation, phases=[water])
