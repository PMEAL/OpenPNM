# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 14:59:50 2019

@author: work
"""
# sample code for relperm
import openpnm as op
import numpy as np

pn = op.network.Cubic(shape=[20, 20, 20], spacing=0.00006)
proj = pn.project
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

##initially fully saturated with water
#########plt.hist(phys_oil['throat.entry_pressure'],bins=100)
##to see the difference
##########plt.hist(phys_oil['throat.entry_pressure'],bins=100)
#inv=op.algorithms.InvasionPercolation(phase=oil,network=pn,project=proj)
#inv.setup(phase=oil,entry_pressure='throat.entry_pressure',pore_volume='pore.volume', throat_volume='throat.volume')
#inlets = pn.pores(['top'])
#used_inlets = [inlets[x] for x in range(0, len(inlets), 2)]
#outlets = pn.pores(['bottom'])
#used_outlets = [outlets[x] for x in range(0, len(outlets), 2)]
#inv.set_inlets(pores=used_inlets)
##run without apply trapping (disconnected clusters grow until touch a connected cluster and stop then)
#inv.run()
##bb=inv.results(Snwp=0.2)
##for invaded zone, find max entry pressure and save it
##have the true index for pore.occupancy and throat.occupancy
##see from their physics model their entry pressure
##get the maximum value of that pressure
##there you go
#"""
#res1=inv.results(Snwp=1)
#occ_ps=res1['pore.occupancy']
#occ_ts=res1['throat.occupancy']
#max_pthroat=np.max(phys_oil['throat.entry_pressure'][occ_ts])
#"""
#Snwparr =  []
#Pcarr =  []
#Sarr=np.linspace(0,1,num=6)
#for Snw in Sarr:
#    res1=inv.results(Snwp=Snw)
#    occ_ts=res1['throat.occupancy']
#    if np.any(occ_ts):
#        max_pthroat=np.max(phys_oil['throat.entry_pressure'][occ_ts])
#        Pcarr.append(max_pthroat)
#        Snwparr.append(Snw)
#        
#y=np.array(Pcarr[:])
#x=1.0-np.array(Snwparr[:])
####plt.xticks(np.arange(x.min(), x.max(), 0.05))
####plt.yticks(np.arange(y.min(), y.max(),0.1))
####plt.plot(x, y)
####plt.xlabel('Invading Phase Saturation')
####plt.ylabel('Capillary Pressure')
####plt.grid(True)  
#  
#filter_pc=y
#sat=1-x
#oil['pore.invasion_sequence']=inv['pore.invasion_sequence']
##proj.export_data(filename='pore_proj',filetype='xdmf')
##################op.io.VTK.save(network=pn, phases=[water,oil])
Hinlets = [pn.pores(['top']), pn.pores(['top']),
              pn.pores(['top'])]
inlets=[]
outlets = [pn.pores(['bottom']), pn.pores(['bottom']),
              pn.pores(['bottom'])]
inv=op.algorithms.InvasionPercolation(phase=oil,network=pn)
inv.setup(phase=oil,entry_pressure='throat.entry_pressure',pore_volume='pore.volume', throat_volume='throat.volume')
pore_inv_seq=[]
throat_inv_seq=[]

def update_phase_and_phys(satnum, invnum, pore_occ, throat_occ):
    oil['pore.occupancy'] = pore_occ[str(invnum)][satnum]
    water['pore.occupancy'] = 1-pore_occ[str(invnum)][satnum]
    oil['throat.occupancy'] = throat_occ[str(invnum)][satnum]
    water['throat.occupancy'] = 1-throat_occ[str(invnum)][satnum]
    print('oil pore', oil['pore.occupancy'])
    print('water pore', water['pore.occupancy'])
    #adding multiphase conductances
    mode='strict'
    water.add_model(model=op.models.physics.multiphase.conduit_conductance,
                       propname='throat.conduit_hydraulic_conductance',
                       throat_conductance='throat.hydraulic_conductance',
                       mode=mode)
    oil.add_model(model=op.models.physics.multiphase.conduit_conductance,
                         propname='throat.conduit_hydraulic_conductance',
                         throat_conductance='throat.hydraulic_conductance',
                         mode=mode)
    
perm_water = {'0': [],'1': [],'2': []}
perm_oil = {'0': [],'1': [],'2': []}
pore_occ = {'0': [],'1': [],'2': []}
throat_occ = {'0': [],'1': [],'2': []}
single_perms_water = [None,None,None]
single_perms_oil = [None,None,None]

[amax, bmax, cmax] = np.max(pn['pore.coords'], axis=0)
[amin, bmin, cmin] = np.min(pn['pore.coords'], axis=0)
lx = amax-amin
ly = bmax-bmin
lz = cmax-cmin
da = lx*ly
dl = lz

def top_b(lx,ly,lz):
    da = lx*ly
    dl = lz
    res_2=[da,dl]
    return res_2

def left_r(lx,ly,lz):
    
    da = lx*lz
    dl = ly
    res_2=[da,dl]
    return res_2

def front_b(lx,ly,lz):
    da = ly*lz
    dl = lx
    res_2=[da,dl]
    return res_2
options = {0 : top_b(lx,ly,lz),1 : left_r(lx,ly,lz),2 : front_b(lx,ly,lz)}

for i in range(len(Hinlets)):
    inlets.append([Hinlets[i][x] for x in range(0, len(Hinlets[i]), 2)])
for i in range(len(inlets)):
    inv.set_inlets(pores=inlets[i])
    inv.run()
    BC1_pores = Hinlets[i]
    BC2_pores = outlets[i]
    [da,dl]=options[i]
    #kw
    Stokes_alg_single_phase_water = op.algorithms.StokesFlow(network=pn, phase=water)
    Stokes_alg_single_phase_water.setup(conductance='throat.hydraulic_conductance')
    Stokes_alg_single_phase_water._set_BC(pores=BC1_pores, bctype='value', bcvalues=100000)
    Stokes_alg_single_phase_water._set_BC(pores=BC2_pores, bctype='value', bcvalues=1000)
#    Stokes_alg_single_phase_water.set_value_BC(values=10**4, pores=BC1_pores)
#    Stokes_alg_single_phase_water.set_value_BC(values=100, pores=BC2_pores)
    Stokes_alg_single_phase_water.run()
    single_perms_water[i] = Stokes_alg_single_phase_water.calc_effective_permeability(domain_area=da, domain_length=dl,inlets=BC1_pores, outlets=BC2_pores)
    proj.purge_object(obj=Stokes_alg_single_phase_water)
    Stokes_alg_single_phase_oil = op.algorithms.StokesFlow(network=pn, phase=oil)
    Stokes_alg_single_phase_oil.setup(conductance='throat.hydraulic_conductance')
    Stokes_alg_single_phase_oil._set_BC(pores=BC1_pores, bctype='value', bcvalues=10000)
    Stokes_alg_single_phase_oil._set_BC(pores=BC2_pores, bctype='value', bcvalues=1000)
    Stokes_alg_single_phase_oil.run()
    single_perms_oil[i] = Stokes_alg_single_phase_oil.calc_effective_permeability(domain_area=da, domain_length=dl,inlets=BC1_pores, outlets=BC2_pores)
    proj.purge_object(obj=Stokes_alg_single_phase_oil)
    pore_inv_seq=inv['pore.invasion_sequence']
    throat_inv_seq=inv['throat.invasion_sequence']
    Snwparr =  []
    Pcarr =  []
    Sarr=np.linspace(0,1,num=15)
    for Snw in Sarr:
        res1=inv.results(Snwp=Snw)
        occ_ts=res1['throat.occupancy']
        if np.any(occ_ts):
            max_pthroat=np.max(phys_oil['throat.entry_pressure'][occ_ts])
            Pcarr.append(max_pthroat)
            Snwparr.append(Snw)
    c=-1
    for Sp in Snwparr:
        c=c+1
        res=inv.results(Sp)
        pore_occ[str(i)].append(res['pore.occupancy'])
        throat_occ[str(i)].append(res['throat.occupancy'])
        update_phase_and_phys(c, i, pore_occ, throat_occ)
        print('sat is equal to', Sp)
        Stokes_alg_multi_phase_water = op.algorithms.StokesFlow(network=pn,phase=water)
        Stokes_alg_multi_phase_water.setup(conductance='throat.conduit_hydraulic_conductance')
        Stokes_alg_multi_phase_water.set_value_BC(values=100000, pores=BC1_pores)
        Stokes_alg_multi_phase_water.set_value_BC(values=1000, pores=BC2_pores)
        #oil
        Stokes_alg_multi_phase_oil = op.algorithms.StokesFlow(network=pn,phase=oil)
        Stokes_alg_multi_phase_oil.setup(conductance='throat.conduit_hydraulic_conductance')
        Stokes_alg_multi_phase_oil.set_value_BC(values=100000, pores=BC1_pores)
        Stokes_alg_multi_phase_oil.set_value_BC(values=1000, pores=BC2_pores)
        # Run Multiphase algs
        Stokes_alg_multi_phase_water.run()
        Stokes_alg_multi_phase_oil.run()
        effective_permeability_water_multi = Stokes_alg_multi_phase_water.calc_effective_permeability(domain_area=da, domain_length=dl)
        effective_permeability_oil_multi = Stokes_alg_multi_phase_oil.calc_effective_permeability(domain_area=da, domain_length=dl)
        relative_eff_perm_water = effective_permeability_water_multi/single_perms_water[i]
        relative_eff_perm_oil = effective_permeability_oil_multi/single_perms_oil[i]
        perm_water[str(i)].append(relative_eff_perm_water)
        perm_oil[str(i)].append(relative_eff_perm_oil)
        proj.purge_object(obj=Stokes_alg_multi_phase_water)
        proj.purge_object(obj=Stokes_alg_multi_phase_oil)





