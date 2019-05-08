# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 14:59:50 2019

@author: work
"""
# sample code for relperm
import openpnm as op
import matplotlib.pyplot as plt
import numpy as np
pn = op.network.Cubic(shape=[20, 20, 20], spacing=0.00006)
geom = op.geometry.StickAndBall(network=pn, pores=pn['pore.all'],
                                throats=pn['throat.all'])
oil = op.phases.GenericPhase(network=pn, name='oil')
water = op.phases.GenericPhase(network=pn, name='water')
oil['pore.viscosity']=0.547
oil['throat.surface_tension'] = 0.072
oil['pore.surface_tension']=0.072
oil['pore.contact_angle']=110
water['throat.contact_angle'] = 70
water['pore.contact_angle'] = 70
water['throat.surface_tension'] = 0.0483
water['pore.surface_tension'] = 0.0483
water['pore.viscosity']=0.4554
mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
oil.add_model(propname='throat.hydraulic_conductance',
              model=mod)
oil.add_model(propname='throat.entry_pressure',
              model=op.models.physics.capillary_pressure.washburn)
water.add_model(propname='throat.hydraulic_conductance',
                model=mod)
water.add_model(propname='throat.entry_pressure',
                model=op.models.physics.capillary_pressure.washburn)
Finlets_init = {'x': pn.pores('left'), 
                'y': pn.pores('front'),
                'z': pn.pores('top')}
Finlets=dict()
for key in Finlets_init.keys():
    Finlets.update({key: ([Finlets_init[key][x] for x in \
                           range(0, len(Finlets_init[key]), 2)])})
Foutlets_init = {'x': pn.pores('right'), 
                   'y': pn.pores('back'),
                   'z': pn.pores('bottom')}
ip = op.algorithms.InvasionPercolation(network=pn, phase=oil)
ip.set_inlets(pores=Finlets['x'])
ip.run()
rp = op.algorithms.DirectionalRelativePermeability(network=pn)
rp.setup(invading_phase=oil, defending_phase=water,
         pore_invasion_sequence=ip['pore.invasion_sequence'],
         throat_invasion_sequence=ip['throat.invasion_sequence'])
rp.set_inlets(pores=Finlets_init['y'])
rp.set_outlets(pores=Foutlets_init['y'])
rp.run(Snw_num=50, IP_pores=Finlets['y'])
results=rp.get_Kr_data()
rp.plot_Kr_curve()
