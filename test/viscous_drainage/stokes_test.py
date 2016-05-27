# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 10:22:44 2016

@author: stadelmanm
"""
import scipy as sp
import OpenPNM
import OpenPNM.Utilities.IO as io
from OpenPNM.Network import tools
#
workspace = OpenPNM.Base.Workspace()
workspace.loglevel = 'WARNING'
#
xdim = 11
ydim = 11
pn = OpenPNM.Network.Cubic(shape=[xdim, ydim, 3], spacing=0.0001,connectivity=18)
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
#
# setting geometery for entire network to ensure proper generation
geom = OpenPNM.Geometry.Toray090(network=pn,
                                 pores=pn.Ps,
                                 throats=pn.Ts)
#
# replacing Toray090 with Boundary class in appropriate locations
geom.set_locations(pores=pn.pores(labels='*_boundary'),
                   mode='remove')
geom.set_locations(throats=pn.throats(labels='*_boundary'),
                   mode='remove')
bound_geom = OpenPNM.Geometry.Boundary(network=pn,
                                       pores=pn.pores(labels='*_boundary'),
                                       throats=pn.throats(labels='*_boundary'))



water = OpenPNM.Phases.Water(network=pn)
phys = OpenPNM.Physics.Standard(network=pn,phase=water,pores=pn.Ps,
                                throats=pn.Ts)
#
Stokes = OpenPNM.Algorithms.StokesFlow(network=pn, phase=water)
Stokes.set_boundary_conditions(bctype='Dirichlet',
                               bcvalue=0.00,
                               pores=pn.pores('back_boundary'))
Stokes.set_boundary_conditions(bctype='Neumann_group',
                               mode='merge',
                               bcvalue=(5.0/3.0)*-1e-11,
                               pores=pn.pores('front_boundary'))
Stokes.setup(phase=water)
Stokes.run()
Stokes.return_results()
#
all_pores_q = Stokes.rate(pores = [0,7], mode='single')
q_in = Stokes.rate(pores = pn.pores('front_boundary'))
q_out = Stokes.rate(pores = pn.pores('back_boundary'))
dp = sp.average(water['pore.pressure'][pn.pores('front_boundary')])
print(q_in,q_out,dp)
#
io.VTK.save(pn,'stokes_pn',[water])



