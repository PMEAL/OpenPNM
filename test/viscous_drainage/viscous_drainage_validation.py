# -*- coding: utf-8 -*-
"""
Created on Mon May 16 16:09:47 2016

@author: stadelmanm
"""
import scipy as sp
import OpenPNM
import OpenPNM.Utilities.IO as io
from OpenPNM.Network import tools
from OpenPNM.Geometry import models as gm
from OpenPNM.Algorithms.__InvasionPercolation__ import InvasionPercolation
from OpenPNM.Algorithms.__ViscousDrainage__ import ViscousDrainage
#
#
workspace = OpenPNM.Base.Workspace()
workspace.loglevel = 'WARNING'
xdim = 5
ydim = 5
pn = OpenPNM.Network.Cubic(shape=[xdim, ydim, 3], spacing=0.0001,connectivity=18,
                           name='vd_validate_net')
#
# trimming top and bottom pores to reduce it down to a 2-D network
tools.trim(pn,pores=pn.pores(labels=['top','bottom']))
#
# trimming every other pore to produce a diamond shapped lattice
pn['pore.to_trim'] = False
for i in range(ydim):
    mod = i % 2
    for j in range(xdim):
        if (j % 2 == mod):
            pn['pore.to_trim'][i*ydim + j] = True
tools.trim(pn,pores=pn.pores(labels=['to_trim']))
#
# adding boundary nodes
pn['throat.internal'] = sp.ones(pn.Nt,dtype=bool)
pn.add_boundaries()
tools.trim(pn,pores=pn.pores(labels=['top_boundary','bottom_boundary',
                                     'left_boundary','right_boundary']))
# setting geometry for entire network to ensure proper generation
geom = OpenPNM.Geometry.Toray090(network=pn,pores=pn.Ps,throats=pn.Ts,
                                 name='toray_geom')
#
geom['pore.seed'] = gm.misc.random(pn.Np,seed=7271993)
geom['throat.seed'] = gm.misc.random(pn.Nt,seed=5181993)
geom.regenerate()
#
# replacing Toray090 with Boundary in appropriate locations
geom.set_locations(pores=pn.pores(labels='*_boundary'),mode='remove')
geom.set_locations(throats=pn.throats(labels='*_boundary'),mode='remove')
bound_geom = OpenPNM.Geometry.Boundary(network=pn,
                                       pores=pn.pores(labels='*_boundary'),
                                       throats=pn.throats(labels='*_boundary'),
                                       name='bound_geom')
#
air = OpenPNM.Phases.Air(network=pn, name='air')
water = OpenPNM.Phases.Water(network=pn, name='water')
#
phys_air = OpenPNM.Physics.Standard(network=pn,phase=air,pores=pn.Ps,
                                    throats=pn.Ts,name='phys_air')
phys_water = OpenPNM.Physics.Standard(network=pn,phase=water,pores=pn.Ps,
                                      throats=pn.Ts,name='phys_water')
model = OpenPNM.Physics.models.capillary_pressure.washburn
phys_water.add_model(propname='throat.capillary_pressure',model=model)
#
# setting bulk injection rate
q = 1.0E-3 / 1.0E6 / 60.0 #convertong 1ml/min to m^3/sec
#
# setting up a stable displacement simulation
VD = ViscousDrainage(network=pn, phase=air)
VD.setup(invading_phase=water,defending_phase=air,injection_rate=q,
         max_steps=1000)
#
VD.set_inlets(pores=pn.pores('front_boundary'))
VD.set_outlets(pores=pn.pores('back_boundary'))
VD.run()
VD.return_results()
#print(VD.rate(pn.pores(),mode='single',phase='invading'))
#print(VD._th_q)
#
pn['pore.number'] = sp.arange(0,pn.Np)
pn['throat.number'] = sp.arange(0,pn.Nt)
io.VTK.save(pn,'temp_files/viscous_drainage_pn',[water,air])
#
workspace.purge_object(pn,mode='complete')
