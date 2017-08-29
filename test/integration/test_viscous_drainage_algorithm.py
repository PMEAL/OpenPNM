import scipy as sp
import OpenPNM
from OpenPNM.Network import tools
from OpenPNM.Geometry import models as gm
from OpenPNM.Algorithms.__ViscousDrainage__ import ViscousDrainage
#
mgr = OpenPNM.Base.Workspace()
mgr.loglevel = 60
#
xdim = 5
ydim = 5
pn = OpenPNM.Network.Cubic(shape=[xdim, ydim, 3],
                           spacing=0.0001,
                           connectivity=18)
#
# trimming top and bottom pores to reduce it down to a 2-D network
tools.trim(pn, pores=pn.pores(labels=['top', 'bottom']))
#
# trimming every other pore to produce a diamond shapped lattice
pn['pore.to_trim'] = False
for i in range(ydim):
    mod = i % 2
    for j in range(xdim):
        if (j % 2 == mod):
            pn['pore.to_trim'][i*ydim + j] = True
tools.trim(pn, pores=pn.pores(labels=['to_trim']))
#
# adding boundary nodes
pn['throat.internal'] = sp.ones(pn.Nt, dtype=bool)
pn.add_boundaries()
tools.trim(pn, pores=pn.pores(labels=['top_boundary', 'bottom_boundary',
                                      'left_boundary', 'right_boundary']))
#
# setting geometry for entire network to ensure proper generation
geom = OpenPNM.Geometry.Toray090(network=pn, pores=pn.Ps, throats=pn.Ts)
#
# regnerating with set seed for testing
geom['pore.seed'] = gm.misc.random(pn.Np, seed=7271993)
geom['throat.seed'] = gm.misc.random(pn.Nt, seed=5181993)
geom.regenerate()
#
# replacing Toray090 with Boundary in appropriate locations
geom.set_locations(pores=pn.pores(labels='*_boundary'), mode='remove')
geom.set_locations(throats=pn.throats(labels='*_boundary'), mode='remove')
bound_geom = OpenPNM.Geometry.Boundary(network=pn,
                                       pores=pn.pores(labels='*_boundary'),
                                       throats=pn.throats(labels='*_boundary'))
#
# creating phases and physics
air = OpenPNM.Phases.Air(network=pn)
water = OpenPNM.Phases.Water(network=pn)
phys_air = OpenPNM.Physics.Standard(network=pn, phase=air, pores=pn.Ps,
                                    throats=pn.Ts)
phys_water = OpenPNM.Physics.Standard(network=pn, phase=water, pores=pn.Ps,
                                      throats=pn.Ts)
model = OpenPNM.Physics.models.capillary_pressure.washburn
phys_water.add_model(propname='throat.capillary_pressure', model=model)
#
# Running algorithm
run_args = {
    'injection_rate': 1.0E-3 / 1.0E6 / 60.0,
    'sat_tol': 1.0E-6,
    'max_steps': 10000,
    'exit_on_breakthough': True
}
VD = ViscousDrainage(network=pn, phase=air)
VD.setup(invading_phase=water, **run_args)
VD.set_inlets(pores=pn.pores('front_boundary'))
VD.set_outlets(pores=pn.pores('back_boundary'))
VD.run()
VD.return_results()
#
#
assert VD.sim_stats['break_through_step'] > 0
