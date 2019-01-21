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
oil = op.phases.GenericPhase(network=pn)
water = op.phases.GenericPhase(network=pn)
oil['pore.viscosity']=0.547
oil['throat.contact_angle'] =110
oil['throat.surface_tension'] = 0.072
oil['pore.surface_tension']=0.072
oil['pore.contact_angle']=110
water['throat.contact_angle'] = 70
water['pore.contact_angle'] = 70
water['throat.surface_tension'] = 0.0483
water['pore.surface_tension'] = 0.0483
water['pore.viscosity']=0.4554
phys_water= op.physics.GenericPhysics(network=pn, phase=water, geometry=geom)
phys_oil = op.physics.GenericPhysics(network=pn, phase=oil, geometry=geom)
mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys_oil.add_model(propname='throat.hydraulic_conductance',
                              model=mod)
phys_oil.add_model(propname='throat.entry_pressure',
                              model=op.models.physics.capillary_pressure.washburn)
phys_water.add_model(propname='throat.hydraulic_conductance',
                              model=mod)
phys_water.add_model(propname='throat.entry_pressure',
                              model=op.models.physics.capillary_pressure.washburn)
inv=op.algorithms.InvasionPercolation(phase=oil,network=pn)
inv.setup(phase=oil,entry_pressure='throat.entry_pressure',pore_volume='pore.volume', throat_volume='throat.volume')
inlets = [pn.pores(['top']), pn.pores(['top']),
              pn.pores(['top'])]
outlets = [pn.pores(['bottom']), pn.pores(['bottom']),
              pn.pores(['bottom'])]
num=15
plt.figure(1)
final={'sat': [], 'K_inv': [], 'K_def': []}
for i in range(len(inlets)):
    inv.set_inlets(pores=inlets[i])
    inv.run()
    Snwparr =  []
    Pcarr =  []
    Sarr=np.linspace(0,1,num=num)
    for Snw in Sarr:
        res1=inv.results(Snwp=Snw)
        occ_ts=res1['throat.occupancy']
        if np.any(occ_ts):
            max_pthroat=np.max(phys_oil['throat.entry_pressure'][occ_ts])
            Pcarr.append(max_pthroat)
            Snwparr.append(Snw)
    pore_occ=[]
    throat_occ=[]
    for Sp in Snwparr:
        res=inv.results(Sp)
        pore_occ.append(res['pore.occupancy'])
        throat_occ.append(res['throat.occupancy'])
    relcalc=op.algorithms.RelativePermeability(network=pn)
    relcalc.setup(inv_phase=oil, def_phase=water, sats=Snwparr,
              pore_inv_seq=inv['pore.invasion_sequence'],
              throat_inv_seq=inv['throat.invasion_sequence'],
              inlets=inlets[i],
              outlets=outlets[i],
              pore_occ=pore_occ,
              throat_occ=throat_occ)
    results=relcalc.run()
    final['sat'].append(results['sat'])
    final['K_inv'].append(results['K_rel_inv'])
    final['K_def'].append(results['K_rel_def'])
