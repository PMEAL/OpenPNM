import OpenPNM
import matplotlib.pyplot as plt

Lc = 40.5e-6

#1 setting up network
sgl = OpenPNM.Network.Cubic([10, 10, 10], spacing=Lc, name='SGL10BA')
sgl.add_boundaries()

#2 set up geometries
Ps = sgl.pores('boundary',mode='difference')
Ts = sgl.find_neighbor_throats(pores=Ps,mode='intersection',flatten=True)
geo = OpenPNM.Geometry.SGL10(network=sgl,pores=Ps,throats=Ts,name='geo')

Ps = sgl.pores('boundary')
Ts = sgl.find_neighbor_throats(pores=Ps,mode='not_intersection')
boun = OpenPNM.Geometry.Boundary(network=sgl,pores=Ps,throats=Ts,name='boun')

#constrict throats in the y and z directions
throats = sgl.throats('geo')
connected_pores = sgl.find_connected_pores(throats)
x1 = [sgl['pore.coords'][pair[0]][0] for pair in connected_pores]
x2 = [sgl['pore.coords'][pair[1]][0] for pair in connected_pores]
same_x = [x - y == 0 for x, y in zip(x1,x2)]
factor = [s*.95 + (not s)*1 for s in same_x]
throat_diameters = sgl['throat.diameter'][throats]*factor
geo['throat.diameter']=throat_diameters

#remove the regeneration ability of the diameter pore and throat properties
geo.remove_model('pore.diameter')
geo.remove_model('throat.diameter')
boun.remove_model('pore.diameter')
boun.remove_model('throat.diameter')

#reset aspects relying on pore and throat sizes
geo.regenerate()
boun.regenerate()

#set up phases
air = OpenPNM.Phases.Air(network = sgl, name = 'air')
water = OpenPNM.Phases.Water(network = sgl, name = 'water')

#calculating all phase values
air.regenerate()
water.regenerate()

#reset pore contact angle
water['pore.contact_angle'] = 100

#1 create physics objects associated with our phases
Ps = sgl.pores()
Ts = sgl.throats()
phys_water = OpenPNM.Physics.Standard(network=sgl,phase=water,pores=Ps,throats=Ts,dynamic_data=True,name='standard_water_physics')
phys_air = OpenPNM.Physics.Standard(network=sgl,phase=air,pores=Ps,throats=Ts,dynamic_data=True,name='standard_air_physics')

#2 calculating physics properties (capillary pressure, hydraulic conductance, etc)
phys_water.regenerate()
phys_air.regenerate()

inlets = sgl.pores('bottom_boundary')
used_inlets = [inlets[x] for x in range(0, len(inlets), 2)]

#using every other pore in the bottom and boundary as an inlet
#prevents extremely small diffusivity and permeability values in the z direction
used_inlets = [inlets[x] for x in range(0, len(inlets), 2)]

OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=sgl,loglevel=30,invading_phase = water, defending_phase = air)
OP_1.run(inlets = used_inlets,npts=100)

sat = []
perm_air = {'0': [], '1': [], '2': []}
diff_air = {'0': [], '1': [], '2': []}
perm_water = {'0': [], '1': [], '2': []}
diff_water = {'0': [], '1': [], '2': []}

max_inv_seq = max(OP_1['throat.inv_seq'])

num_seq = 20
for x in range(num_seq+1):
    OP_1.return_results(sat = x/num_seq)

    #printing out so we know how far along we are
    print('seq = '+str(round(max_inv_seq*(x/num_seq)))+' Seq out of '+str(round(max_inv_seq))+' total sequences')

    final_pores = water['pore.occupancy']
    pore_volumes = sgl['pore.volume']
    final_throats = water['throat.occupancy']
    throat_volumes = sgl['throat.volume']

    saturation = (sum(final_pores*pore_volumes) + sum(final_throats*throat_volumes))/(sum(pore_volumes) + sum(throat_volumes))

    sat.append(saturation)

    #adding multiphase conductances
    phys_air.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
               propname='throat.conduit_diffusive_conductance',
               throat_conductance='throat.diffusive_conductance')
    phys_water.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
               propname='throat.conduit_diffusive_conductance',
               throat_conductance='throat.diffusive_conductance')
    phys_air.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
               propname='throat.conduit_hydraulic_conductance',
               throat_conductance='throat.hydraulic_conductance')
    phys_water.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
               propname='throat.conduit_hydraulic_conductance',
               throat_conductance='throat.hydraulic_conductance')

    bounds = [['front', 'back'], ['left', 'right'], ['top', 'bottom']]
    
    for bound_increment in range(len(bounds)):

        #run Stokes Flow and find Permeability
        #single phase
        Stokes_alg_single_phase_air = OpenPNM.Algorithms.StokesFlow(name='Stokes_alg_single_phase_air'+str(bound_increment)+str(x),network=sgl,phase=air)
        Stokes_alg_single_phase_water = OpenPNM.Algorithms.StokesFlow(name='Stokes_alg_single_phase_water'+str(bound_increment)+str(x),network=sgl,phase=water)
        
        Fickian_alg_single_phase_air = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_single_phase_air'+str(bound_increment)+str(x),network=sgl,phase=air)
        Fickian_alg_single_phase_water = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_single_phase_water'+str(bound_increment)+str(x),network=sgl,phase=water)
        
        Stokes_alg_multi_phase_air = OpenPNM.Algorithms.StokesFlow(name='Stokes_alg_multi_phase_air'+str(bound_increment)+str(x),network=sgl,phase=air)
        Stokes_alg_multi_phase_water = OpenPNM.Algorithms.StokesFlow(name='Stokes_alg_multi_phase_water'+str(bound_increment)+str(x),network=sgl,phase=water)
        
        Fickian_alg_multi_phase_air = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_multi_phase_air'+str(bound_increment)+str(x),network=sgl,phase=air)
        Fickian_alg_multi_phase_water = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_multi_phase_water'+str(bound_increment)+str(x),network=sgl,phase=water)
        
        BC1_pores = sgl.pores(labels=bounds[bound_increment][0]+'_boundary')
        BC2_pores = sgl.pores(labels=bounds[bound_increment][1]+'_boundary')

        #BC1
        Stokes_alg_single_phase_air.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.6,pores=BC1_pores)
        Stokes_alg_single_phase_water.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.6,pores=BC1_pores)
        Fickian_alg_single_phase_air.set_boundary_conditions(bctype='Dirichlet',bcvalue=.6,pores=BC1_pores)
        Fickian_alg_single_phase_water.set_boundary_conditions(bctype='Dirichlet',bcvalue=.6,pores=BC1_pores)
        
        Stokes_alg_multi_phase_air.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.6,pores=BC1_pores)
        Stokes_alg_multi_phase_water.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.6,pores=BC1_pores)
        Fickian_alg_multi_phase_air.set_boundary_conditions(bctype='Dirichlet',bcvalue=.6,pores=BC1_pores)
        Fickian_alg_multi_phase_water.set_boundary_conditions(bctype='Dirichlet',bcvalue=.6,pores=BC1_pores)
        
        #BC2
        Stokes_alg_single_phase_air.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.2,pores=BC2_pores)
        Stokes_alg_single_phase_water.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.2,pores=BC2_pores)
        Fickian_alg_single_phase_air.set_boundary_conditions(bctype='Dirichlet',bcvalue=.2,pores=BC2_pores)
        Fickian_alg_single_phase_water.set_boundary_conditions(bctype='Dirichlet',bcvalue=.2,pores=BC2_pores)
        
        Stokes_alg_multi_phase_air.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.2,pores=BC2_pores)
        Stokes_alg_multi_phase_water.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.2,pores=BC2_pores)
        Fickian_alg_multi_phase_air.set_boundary_conditions(bctype='Dirichlet',bcvalue=.2,pores=BC2_pores)
        Fickian_alg_multi_phase_water.set_boundary_conditions(bctype='Dirichlet',bcvalue=.2,pores=BC2_pores)
        
        #run algorithms with proper conduit conductance
        Stokes_alg_single_phase_air.run(conductance = 'hydraulic_conductance')
        Stokes_alg_single_phase_water.run(conductance = 'hydraulic_conductance')
        Fickian_alg_single_phase_air.run(conductance = 'diffusive_conductance')
        Fickian_alg_single_phase_water.run(conductance = 'diffusive_conductance')
        
        Stokes_alg_multi_phase_air.run(conductance = 'conduit_hydraulic_conductance')
        Stokes_alg_multi_phase_water.run(conductance = 'conduit_hydraulic_conductance')
        Fickian_alg_multi_phase_air.run(conductance = 'conduit_diffusive_conductance')
        Fickian_alg_multi_phase_water.run(conductance = 'conduit_diffusive_conductance')
        
        #calc effective properties
        effective_permeability_air_single = Stokes_alg_single_phase_air.calc_eff_permeability()  
        effective_diffusivity_air_single = Fickian_alg_single_phase_air.calc_eff_diffusivity()
        effective_permeability_water_single = Stokes_alg_single_phase_water.calc_eff_permeability()  
        effective_diffusivity_water_single = Fickian_alg_single_phase_water.calc_eff_diffusivity()
        
        effective_permeability_air_multi = Stokes_alg_multi_phase_air.calc_eff_permeability()  
        effective_diffusivity_air_multi = Fickian_alg_multi_phase_air.calc_eff_diffusivity()
        effective_permeability_water_multi = Stokes_alg_multi_phase_water.calc_eff_permeability()  
        effective_diffusivity_water_multi = Fickian_alg_multi_phase_water.calc_eff_diffusivity()
        
        relative_eff_perm_air = effective_permeability_air_multi/effective_permeability_air_single
        relative_eff_perm_water = effective_permeability_water_multi/effective_permeability_water_single
        relative_eff_diff_air = effective_diffusivity_air_multi/effective_diffusivity_air_single
        relative_eff_diff_water = effective_diffusivity_water_multi/effective_diffusivity_water_single

        perm_air[str(bound_increment)].append(relative_eff_perm_air)
        diff_air[str(bound_increment)].append(relative_eff_diff_air)
        perm_water[str(bound_increment)].append(relative_eff_perm_water)
        diff_water[str(bound_increment)].append(relative_eff_diff_water)



from matplotlib.font_manager import FontProperties

#Data points taken directly from Gostick's graphs using GraphClick
gostick_saturation_1 = [0.008, 0.04, 0.093, 0.14, 0.193, 0.246, 0.293, 0.337, 0.395, 0.442, 0.496,
                        0.542, 0.59, 0.641, 0.687, 0.748, 0.793, 0.838, 0.894, 0.945, 0.986]
gostick_perm_air_case1 = [0.917, 0.821, 0.68, 0.568, 0.466, 0.366, 0.286, 0.204, 0.144, 0.096, 0.051, 0.024,
                          0.003, -1.08E-04, -1.96E-04, -3.12E-04, -3.97E-04, -4.84E-04, -5.90E-04, 0.002, 0.002]
gostick_saturation_2 = [0.99, 0.899, 0.847, 0.802, 0.75, 0.701, 0.645, 0.594, 0.546, 0.497, 0.449,
                        0.398, 0.348, 0.298, 0.245, 0.196, 0.147, 0.094, 0.044, 0.003]
gostick_perm_water = [0.935, 0.774, 0.709, 0.664, 0.618, 0.572, 0.514, 0.461, 0.401, 0.347,
                        0.284, 0.211, 0.145, 0.084, 0.044, 0.024, 0.012, 0.001, 0.001, 0.001]

gostick_saturation_3 =[0.006, 0.05, 0.102, 0.151, 0.199, 0.247, 0.297, 0.348, 0.399, 0.447, 0.496,
                    0.546, 0.597, 0.645, 0.699, 0.75, 0.798, 0.846, 0.899, 0.949, 0.983]
gostick_diff_air_case1 = [0.939, 0.836, 0.725, 0.626, 0.531, 0.442, 0.353, 0.27, 0.203, 0.14, 0.085, 0.048,
                          0.008, 5.49E-04, 4.48E-04, 3.50E-04, 2.59E-04, 1.67E-04, 0.003, 0.003, 0.003]
gostick_saturation_4 = [0.985, 0.946, 0.898, 0.846, 0.795, 0.749, 0.695, 0.643, 0.596, 0.545, 0.496, 0.448,
                        0.396, 0.346, 0.298, 0.251, 0.196, 0.146, 0.094]
gostick_diff_water = [0.941, 0.901, 0.853, 0.809, 0.756, 0.7, 0.638, 0.569, 0.503, 0.428, 0.36, 0.291, 0.214, 1.48E-01,
                      8.00E-02, 4.50E-02, 2.30E-02, 1.60E-02, 0.005]

fontP = FontProperties()
fontP.set_size('small')
#setting up subplots
fig = plt.figure(num=1, figsize=(6, 10), dpi=80, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(211)   #top 
ax2 = fig.add_subplot(212)   #bottom 

x_values1 = [x/20 for x in range(21)]
z = '.75'


#plots for subplot1 - strict permeability
p1, = ax1.plot(sat, perm_water['0'], color = 'k', linestyle = '-', marker = 'o')
p2, = ax1.plot(sat, perm_water['1'], color = z, linestyle = '-', marker = 'o')
p3, = ax1.plot(sat, perm_water['2'], color = 'w', linestyle = '-', marker = 'o')
p4, = ax1.plot(sat, perm_air['0'], color = 'k', linestyle = '-', marker = '^')
p5, = ax1.plot(sat, perm_air['1'], color = z, linestyle = '-', marker = '^')
p6, = ax1.plot(sat, perm_air['2'], color = 'w', linestyle = '-', marker = '^')
p10, = ax1.plot(x_values1, [x**(3) for x in x_values1], 'k--')
ax1.plot(x_values1, [(1-x)**(3) for x in x_values1], 'k--')
gs1, = ax1.plot(gostick_saturation_1, gostick_perm_air_case1, color = 'r', linestyle = '-', marker = 'D')
gs2, = ax1.plot(gostick_saturation_2, gostick_perm_water, color = 'r', linestyle = '-', marker = 'o')
ax1.set_ylabel('permeability')
ax1.set_xlabel("saturation")
ax1.set_ylim([0,1])
ax1.set_xlim([0,1])

#need to work on legend to match up with the right things
lgd1 = ax1.legend([p1, p2, p3, p4, p5, p6, p10, gs1, gs2],
           ["KrWater,x", "KrWater,y", "KrWater,z",
           "KrAir,x","KrAir,y","KrAir,z", "a = 3", "Gostick et al \n KrAir,x (case 1)", "Gostick et al \n KrWater,x"], loc='center left', bbox_to_anchor=(1, 0.5), prop = fontP)

#plots for subplot4 - diffusivity
p11, = ax2.plot(sat, diff_water['0'], color = 'k', linestyle = '-', marker = 'o')
p12, = ax2.plot(sat, diff_water['1'], color = z, linestyle = '-', marker = 'o')
p13, = ax2.plot(sat, diff_water['2'], color = 'w', linestyle = '-', marker = 'o')
p14, = ax2.plot(sat, diff_air['0'], color = 'k', linestyle = '-', marker = '^')
p15, = ax2.plot(sat, diff_air['1'], color = z, linestyle = '-', marker = '^')
p16, = ax2.plot(sat, diff_air['2'], color = 'w', linestyle = '-', marker = '^')
p20, = ax2.plot(x_values1, [x**(2) for x in x_values1], 'k--')
ax2.plot(x_values1, [(1-x)**(2) for x in x_values1], 'k--')
gs3, = ax2.plot(gostick_saturation_3, gostick_diff_air_case1, color = 'r', linestyle = '-', marker = 'D')
gs4, = ax2.plot(gostick_saturation_4, gostick_diff_water, color = 'r', linestyle = '-', marker = 'o')
ax2.set_ylabel('diffusivity')
ax2.set_xlabel("saturation")
ax2.set_ylim([0,1])
ax2.set_xlim([0,1])

lgd2 = ax2.legend([p11, p12, p13, p14, p15, p16, p20, gs3, gs4],
           ["DrWater,x", "DrWater,y", "DrWater,z",
           "DrAir,x","DrAir,y","DrAir,z", "a = 2", "Gostick et al \n DrAir,x (case 1)", "Gostick et al \n DrWater,x"], loc='center left', bbox_to_anchor=(1, 0.5), prop = fontP)

fig.subplots_adjust(left=0.13, right=.7, top=0.95, bottom=0.05)

fig.show()
