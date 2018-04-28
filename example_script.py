import openpnm as op
import scipy as sp
ws = op.core.Workspace()

sp.random.seed(0)
pn = op.network.Cubic(shape=[4, 4, 1], spacing=0.0001, name='pn11')
pn.add_boundary_pores()
proj = pn.project

Ps = pn.pores('internal')
Ts = pn.throats('internal')
geom = op.geometry.StickAndBall(network=pn, pores=Ps, throats=Ts,
                                name='intern')

Ps = pn.pores('*boundary')
Ts = pn.throats('*boundary')
boun = op.geometry.StickAndBall(network=pn, pores=Ps, throats=Ts, name='bound')

pn['pore.inlets'] = pn['pore.top_boundary'].copy()
pn['pore.outlets'] = pn['pore.bottom_boundary'].copy()

air = op.phases.Air(network=pn)
water = op.phases.Water(network=pn)
water['throat.viscosity'] = water['pore.viscosity'][0]

mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys_water_1 = op.physics.GenericPhysics(network=pn, phase=water,
                                         geometry=geom)
phys_water_1.add_model(propname='throat.conductance',
                       model=mod,
                       viscosity='throat.viscosity')

phys_water_2 = op.physics.GenericPhysics(network=pn, phase=water,
                                         geometry=boun)
phys_water_2.models = phys_water_1.models.copy()

water['pore.A'] = 1e-10
water['pore.k'] = 2
water.add_model(propname='pore.reaction',
                model=op.models.physics.generic_source_term.standard_kinetics,
                prefactor='pore.A', exponent='pore.k',
                quantity='pore.pressure', regen_mode='normal')

s = {'conductance': 'throat.conductance',
     'quantity': 'pore.pressure'}
alg1 = op.algorithms.GenericTransport(network=pn)
# You can specify phase and settings using the setup method
alg1.setup(phase=water, **s)
alg1.set_dirichlet_BC(pores=pn.pores('inlets'), values=1)
alg1.set_dirichlet_BC(pores=pn.pores('outlets'), values=0)
alg1.run()

# You can also specify phase and settings during initialization
alg2 = op.algorithms.ReactiveTransport(network=pn, phase=water, settings=s)
alg2.settings.update({'tolerance': 0.001})
alg2.set_dirichlet_BC(pores=pn.pores('inlets'), values=1)
alg2.set_source(propname='pore.reaction', pores=pn.pores('outlets'))
alg2.run()
water.update(alg2.results())

alg3 = op.algorithms.TransientReactiveTransport(network=pn, phase=water)
alg3.settings.update(s)
alg3.settings.update({'t_scheme': 'implicit', 'r_tolerance': 1e-6,
                      't_tolerance': 1e-10})
alg3.set_IC(0)
alg3.set_dirichlet_BC(pores=pn.pores('inlets'), values=1)
alg3.set_source(propname='pore.reaction', pores=pn.pores('outlets'))
alg3.run()

#test = alg3['pore.pressure_steady'] - alg2['pore.pressure']
