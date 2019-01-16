# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 14:59:50 2019

@author: work
"""
# sample code for relperm

import openpnm as op
import matplotlib.pyplot as plt
from openpnm import topotools as tt
pn = op.network.Cubic(shape=[3, 3, 3], spacing=0.00006)
proj = pn.project
Psb = pn.pores(['top', 'bottom'])
Tsb = pn.find_neighbor_throats(pores=Psb)
Ps=pn.pores(['top', 'bottom'], mode='not')
Ts = pn.find_neighbor_throats(pores=pn.pores(['top', 'bottom'],
                                             mode='not'), mode='xnor')
a = pn.check_network_health()
tt.trim(network=pn, pores=a['trim_pores'])
geom = op.geometry.StickAndBall(network=pn, pores=pn['pore.all'],
                                throats=pn['throat.all'])
oil = op.phases.GenericPhase(network=pn)
water = op.phases.Water(network=pn)
oil['pore.viscosity']=0.547
oil['throat.contact_angle'] =110
oil['throat.surface_tension'] = 0.072
oil['pore.surface_tension']=0.072
oil['pore.contact_angle']=110
phys_water= op.physics.GenericPhysics(network=pn, phase=water, geometry=geom)
phys_oil = op.physics.GenericPhysics(network=pn, phase=oil, geometry=geom)

mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys_oil.add_model(propname='throat.hydraulic_conductance',
                              model=mod)
phys_oil.add_model(propname='throat.entry_pressure',
                              model=op.models.physics.capillary_pressure.washburn)
# phys_water['pore.entry_pressure'] = -100000
# phys_oil['pore.entry_pressure'] = -1000000
# phys_oil.add_model(propname='pore.entry_pressure',
#                              model=op.models.physics.capillary_pressure.washburn)
phys_water.add_model(propname='throat.hydraulic_conductance',
                              model=mod)
phys_water.add_model(propname='throat.entry_pressure',
                              model=op.models.physics.capillary_pressure.washburn)
relcalc=op.algorithms.RelativePermeability(network=pn)

relcalc.setup(inv_phase=oil, def_phase=water, points=10)
results=relcalc.run()
# aa=relcalc['pore_occ']
# Results = {'k_inv': [], 'k_def': [], 'K_rel_inv': [], 'K_rel_def': []}
x=results['sat']
plt.figure(1)
 for i in range(len(results['k_inv'])):
    y1=results['K_rel_inv'][i]
    y2=results['K_rel_def'][i]
    plt.plot(x,y1)
    plt.plot(x,y2)
#        x=1.0-np.array(Snwparr[:])
#        plt.xticks(np.arange(x.min(), x.max(), 0.05))
#        plt.yticks(np.arange(y.min(), y.max(),0.1))
#        plt.plot(x, y)
#        plt.xlabel('Invading Phase Saturation')
#        plt.ylabel('relative permeability')
