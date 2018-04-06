import openpnm as op
import scipy as sp
ws = op.core.Workspace()

sp.random.seed(0)
pn = op.network.Cubic(shape=[5, 5, 5], spacing=0.0001, name='pn11')
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
                quantity='pore.pressure', regen_mode='deferred')

s = {'conductance': 'throat.conductance',
     'quantity': 'pore.pressure'}
alg1 = op.algorithms.GenericTransport(network=pn, phase=water, settings=s)
alg1.set_dirichlet_BC(pores=pn.pores('inlets'), values=1)
alg1.set_dirichlet_BC(pores=pn.pores('outlets'), values=0)
alg1.run()

alg2 = op.algorithms.ReactiveTransport(network=pn, phase=water, settings=s)
alg2.set_dirichlet_BC(pores=pn.pores('inlets'), values=1)
alg2.set_source_term(propname='pore.reaction', pores=pn.pores('outlets'))
alg2.run()
water.update(alg2.results())

alg3 = op.algorithms.TransientTransport(network=pn, phase=water)
# You can also set the settings afterwards.  Note that some of these
# will have defaults when finally subclassed (i.e. quantity = pressure)
alg3.settings.update({'t_initial': 0,
                      't_final': 1,
                      't_step': 0.25,
                      'conductance': 'throat.conductance',
                      'quantity': 'pore.pressure'})
alg3.set_dirichlet_BC(pores=pn.pores('inlets'), values=1)
alg3.set_IC(values=0)
alg3.run()


alg4 = op.algorithms.TransientReactiveTransport(network=pn, phase=water)
alg4.settings.update(alg3.settings)  # Just copy settings from another alg
alg4.set_dirichlet_BC(pores=pn.pores('inlets'), values=1)
alg4.set_IC(values=0)
alg4.set_source_term(propname='pore.reaction', pores=pn.pores('bottom'))
alg4.run()
