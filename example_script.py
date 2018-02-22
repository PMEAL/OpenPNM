import openpnm as op
import scipy as sp
ws = op.core.Workspace()

sp.random.seed(0)
pn = op.network.Cubic(shape=[5, 5, 5], spacing=0.0001, name='pn11')
pn.add_boundary_pores()
proj = pn.project

geom = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts,
                                settings={'test': 1})

air = op.phases.Air(network=pn)
water = op.phases.Water(network=pn)
water['throat.viscosity'] = water['pore.viscosity'][0]

mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys_water = op.physics.GenericPhysics(network=pn, phase=water, geometry=geom)
phys_water.add_model(propname='throat.conductance',
                     model=mod, viscosity='throat.viscosity')

water['pore.A'] = 1e-10
water['pore.k'] = 2
water.add_model(propname='pore.reaction',
                model=op.models.physics.generic_source_term.standard_kinetics,
                prefactor='pore.A', exponent='pore.k',
                quantity='pore.pressure', regen_mode='deferred')

s = {'conductance': 'throat.conductance',
     'quantity': 'pore.pressure'}
alg1 = op.algorithms.GenericTransport(network=pn, phase=water, settings=s)
alg1.set_dirchlet_BC(pores=pn.pores('top'), values=1)
alg1.set_dirchlet_BC(pores=pn.pores('bottom'), values=0)
alg1.run()

alg2 = op.algorithms.ReactiveTransport(network=pn, phase=water, settings=s)
alg2.set_dirchlet_BC(pores=pn.pores('top'), values=1)
alg2.set_source_term(propname='pore.reaction', pores=pn.pores('bottom'))
alg2.run()
water.update(alg2.results())

alg3 = op.algorithms.TransientTransport(network=pn, phase=water)
# You can also set the settings afterwards.  Note that some of these
# will have defaults when finally subclassed (i.e. quantity = pressure)
alg3.settings.update({'t_initial': 0,
                      't_final': 1,
                      't_step': 0.5,
                      'conductance': 'throat.conductance',
                      'quantity': 'pore.pressure'})
alg3.set_dirchlet_BC(pores=pn.pores('top'), values=1)
alg3.set_IC(values=0)
alg3.run()
