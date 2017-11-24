import openpnm as op
import scipy as sp
import openpnm.geometry.models as gm
ws = op.core.Workspace()

pn = op.network.Cubic(shape=[15, 15, 15], spacing=0.0001)
Ps = pn.pores(['top', 'bottom', 'left', 'right', 'front', 'back'])
pn['pore.surface'] = pn.tomask(pores=Ps)
Ts = pn.find_neighbor_throats(pores=Ps, mode='intersection')
pn['throat.surface'] = pn.tomask(throats=Ts)
geom1 = op.geometry.StickAndBall(network=pn, pores=pn.pores('surface'),
                                 throats=pn.throats('surface'))
geom1['pore.num'] = 1
Ps = pn.pores('surface', mode='not')
pn['pore.internal'] = pn.tomask(pores=Ps)
Ts = pn.throats('surface', mode='not')
pn['throat.internal'] = pn.tomask(throats=Ts)
geom2 = op.geometry.StickAndBall(network=pn, pores=Ps, throats=Ts)
geom2['pore.num'] = 2

air1 = op.phases.Air(network=pn)
air2 = op.phases.Air(network=pn)
mercury = op.phases.Mercury(network=pn, name='Hg')
water1 = op.phases.Water(network=pn)
phys1 = op.physics.GenericPhysics(network=pn, geometry=geom1, phase=mercury)
phys2 = op.physics.GenericPhysics(network=pn, geometry=geom2, phase=mercury)
mercury.add_model(propname='throat.capillary_pressure',
                  model=op.physics.models.capillary_pressure.washburn)
mercury.regenerate_models()

a = pn.simulation

mip = op.algorithms.MIP(network=pn)


pn2 = op.network.Cubic(shape=[6, 6, 6])
pn2.add_model(propname='pore.seed',
              model=gm.pore_misc.random,
              seed=1, num_range=[0, 0.1],
              regen_mode='deferred')
pn2.add_model(propname='pore.diameter',
              model=gm.pore_size.weibull,
              shape=5, scale=0.5, loc=0.1,
              regen_mode='deferred')
pn2.add_model(propname='pore.volume',
              model=gm.pore_volume.sphere,
              regen_mode='deferred')
water2 = op.phases.Water(network=pn2)
water2.add_model(propname='throat.capillary_pressure',
                 model=op.physics.models.capillary_pressure.washburn)

b = pn2.simulation


