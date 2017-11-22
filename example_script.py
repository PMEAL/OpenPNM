import openpnm as op
import scipy as sp
import openpnm.geometry.models as gm
ws = op.core.Workspace()

pn = op.network.Cubic(shape=[5, 5, 5], spacing=1, name='net1')
Ps = pn.pores(['top', 'bottom', 'left', 'right', 'front', 'back'])
pn['pore.surface'] = pn.tomask(pores=Ps)
Ts = pn.find_neighbor_throats(pores=Ps, mode='intersection')
pn['throat.surface'] = pn.tomask(throats=Ts)
geom1 = op.geometry.GenericGeometry(network=pn, name='geo1',
                                    pores=pn.pores('surface'),
                                    throats=pn.throats('surface'))
geom1['pore.num'] = 1
Ps = pn.pores('surface', mode='not')
pn['pore.internal'] = pn.tomask(pores=Ps)
Ts = pn.throats('surface', mode='not')
pn['throat.internal'] = pn.tomask(throats=Ts)
geom2 = op.geometry.StickAndBall(network=pn, name='geo2', pores=Ps, throats=Ts)
geom2['pore.num'] = 2

air1 = op.phases.Air(network=pn, name='air1')
air2 = op.phases.Air(network=pn, name='air2')
mercury = op.phases.Mercury(network=pn, name='Hg')
water1 = op.phases.Water(network=pn, name='water')
phys1 = op.physics.GenericPhysics(network=pn, geometry=geom1, phase=mercury)
phys2 = op.physics.GenericPhysics(network=pn, geometry=geom2, phase=mercury)
phys1['throat.capillary_pressure'] = sp.rand(geom1.Nt)
phys2['throat.capillary_pressure'] = sp.rand(geom2.Nt)

a = pn.simulation

mip = op.algorithms.Drainage(network=pn)
mip.setup(invading_phase=mercury)
mip.set_inlets(pores=pn.pores('surface'))


class CubicWithModels(op.network.Cubic, op.core.ModelsMixin):
    pass


pn2 = CubicWithModels(shape=[6, 6, 6], name='net2')
pn2.add_model(propname='pore.seed',
              model=gm.pore_misc.random,
              seed=1, num_range=[0, 0.1],
              regen_mode='deferred')
pn2.add_model(propname='pore.diameter',
              model=gm.pore_diameter.weibull,
              shape=5, scale=0.5, loc=0.1,
              regen_mode='deferred')
pn2.add_model(propname='pore.volume',
              model=gm.pore_volume.sphere,
              regen_mode='deferred')
water2 = op.phases.Water(network=pn2)
water2.add_model(propname='throat.capillary_pressure',
                 model=op.physics.models.capillary_pressure.washburn)

b = pn2.simulation
