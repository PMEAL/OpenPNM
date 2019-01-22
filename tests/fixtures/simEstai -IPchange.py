import openpnm as op
import numpy as np
import matplotlib.pyplot as plt
import math
import random
mgr = op.Workspace()
mgr.clear()
mgr.keys()
mgr.load_workspace('estail1000.pnm')
pn=mgr['sim_08']['net_01']
proj = pn.project
print(pn)
coord_num_avg=np.mean(pn.num_neighbors(pores=pn.Ps))
del pn['pore.area']
del pn['throat.conduit_lengths.pore1']
del pn['throat.conduit_lengths.pore2']
del pn['throat.conduit_lengths.throat']
del pn['throat.endpoints.tail']
del pn['throat.endpoints.head']
del pn['throat.volume']

pn['pore.diameter']=pn['pore.equivalent_diameter']
pn['throat.diameter']=pn['throat.equivalent_diameter']
"""option for different theta
"""
#"""
Num_pores=len(pn['pore.all'])
Npr=math.floor(Num_pores*3/4)
Num_throats=len(pn['throat.all'])
rndi=random.sample(range(0,Num_pores-1 ),Npr)
WW_pores=np.unique(pn.pores(['all'])[rndi])
WW_throats=pn.find_neighbor_throats(pores=WW_pores)
pn['pore.WW'] = False
pn['pore.WW'][WW_pores]=True
pn['throat.WW']=False
pn['throat.WW'][WW_throats]=True
OW_pores=pn.pores(['pore.WW'],mode='not')
OW_throats=pn.throats(['throat.WW'],mode='not')
pn['pore.OW'] = False
pn['pore.OW'][OW_pores]=True
pn['throat.OW']=False
pn['throat.OW'][OW_throats]=True
#"""
"""end of option for different theta
"""
geom = op.geometry.GenericGeometry(network=pn, pores=pn['pore.all'], throats=pn['throat.all'],project=proj)
#geom=mgr['sim_08']['geo_01']
geom.add_model(propname='throat.endpoints',
                model=op.models.geometry.throat_endpoints.spherical_pores)
geom.add_model(propname='pore.area',
                model=op.models.geometry.pore_area.sphere)
geom.add_model(propname='throat.volume',
                model=op.models.geometry.throat_volume.cylinder)
geom.add_model(propname='throat.conduit_lengths',
                model=op.models.geometry.throat_length.conduit_lengths)
#pn['pore.area']=geom['pore.area']
###plt.hist(pn['pore.diameter'],bins=1000)
###plt.hist(pn['throat.diameter'],bins=1000)
oil = op.phases.GenericPhase(network=pn,project=proj)
water = op.phases.GenericPhase(network=pn,project=proj)
oil['pore.viscosity']=0.547
oil['throat.contact_angle'] =180
oil['throat.surface_tension'] = 0.072
oil['pore.surface_tension']=0.072
oil['pore.contact_angle']=180
water['throat.contact_angle'] = 0
water['pore.contact_angle'] = 0
water['throat.surface_tension'] = 0.0483
water['pore.surface_tension'] = 0.0483
water['pore.viscosity']=0.4554

"""option for different theta
"""
#"""
size=np.random.uniform(0, 60, size=len(pn.pores(['WW'])))

water['pore.contact_angle'][WW_pores]=np.random.uniform(0, 60, size=len(pn.pores(['WW'])))
water['throat.contact_angle'][WW_pores]=np.random.uniform(0, 60, size=len(pn.pores(['WW'])))
water['pore.contact_angle'][OW_pores]=np.random.uniform(110, 160, size=len(pn.pores(['OW'])))
water['throat.contact_angle'][OW_pores]=np.random.uniform(110, 160, size=len(pn.pores(['OW'])))

oil['pore.contact_angle'][WW_pores]=180-water['pore.contact_angle'][WW_pores]
oil['throat.contact_angle'][WW_pores]=180-water['throat.contact_angle'][WW_pores]
oil['pore.contact_angle'][OW_pores]=180-water['pore.contact_angle'][OW_pores]
oil['throat.contact_angle'][OW_pores]=180-water['throat.contact_angle'][OW_pores]
#"""
"""
oil['pore.contact_angle'][WW_pores]=np.random.uniform(110, 160, size=len(pn.pores(['WW'])))
oil['throat.contact_angle'][WW_pores]=np.random.uniform(110, 160, size=len(pn.pores(['WW'])))
oil['pore.contact_angle'][OW_pores]=np.random.uniform(0, 60, size=len(pn.pores(['OW'])))
oil['throat.contact_angle'][OW_pores]=np.random.uniform(0, 60, size=len(pn.pores(['OW'])))
"""
"""end of option for different theta
"""
phys_water= op.physics.GenericPhysics(network=pn, phase=water, geometry=geom,project=proj)
phys_oil = op.physics.GenericPhysics(network=pn, phase=oil, geometry=geom,project=proj)

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

def update_phase_and_phys(satnum, pore_occ, throat_occ):
    oil['pore.occupancy'] = pore_occ[satnum]
    water['pore.occupancy'] = 1-pore_occ[satnum]
    oil['throat.occupancy'] = throat_occ[satnum]
    water['throat.occupancy'] = 1-throat_occ[satnum]
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
    pore_occ=[]
    throat_occ=[]
    c=-1
    for Sp in Snwparr:
        c=c+1
        res=inv.results(Sp)
        pore_occ.append(res['pore.occupancy'])
        throat_occ.append(res['throat.occupancy'])
        update_phase_and_phys(c, pore_occ, throat_occ)
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
#plotting curves
#proj.export_data(filename='pore_proj',filetype='xdmf')
#op.io.VTK.save(network=pn, phases=[water,oil])
#fig = plt.figure(figsize=(8,12), dpi=100, facecolor='w', edgecolor='k')
#sat=(1-sat)
#p1=plt.plot(sat, perm_oil['2'], color = 'r', linestyle = '-', marker = 'o')
#p2=plt.plot(sat, perm_water['2'], color = 'b', linestyle = '-', marker = '^')
#plt.ylabel('Relative Permeability')
#plt.xlabel('SaturationW')
#plt.ylim([0,1])
#plt.xlim([0,1])
##need to work on legend to match up with the right things
#lgd1 = fig.legend([p2, p1],["KrWater", "Kroil"])
#plt.show()






