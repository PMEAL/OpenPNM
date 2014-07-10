# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 16:38:29 2014

@author: Jackie
"""
import OpenPNM
import matplotlib.pyplot as plt


r"""
returns array
array[0] = saturation
array[1] = relative permeability
array[2] = relative effective diffusivity
"""  
Lc = 40.5e-6


#setting up network
=======
sgl = OpenPNM.Network.Cubic(name = 'SGL10BA', loglevel = 40)
sgl.generate(divisions = [26, 26, 10], add_boundaries = True, lattice_spacing = [Lc])

#set up geometries
geo = OpenPNM.Geometry.SGL10(name = 'geo', network = sgl)
geo.set_locations(pores=sgl.pores('internal'),throats='all')

boun = sgl.add_geometry(subclass='Boundary',name='boun')
boun.set_locations(pores=sgl.pores('boundary'))


sgl.regenerate_geometries()

#account for pores that are too big
value = [min(sgl.get_pore_data(prop = 'diameter', locations = x), Lc) for x in geo.pores()]
sgl.set_data(prop='diameter',pores=geo.pores(),data=value)
#account for throats that are too big
value = [min(sgl.get_throat_data(prop = 'diameter', locations = x), Lc) for x in geo.throats()]
sgl.set_data(prop='diameter',throats=geo.throats(),data=value)

#spatial correlation of pore sizes?!

#constricting sgl by .95 in both the z and y direction  
throats = geo.throats()
connected_pores = sgl.find_connected_pores(throats)
x1 = [sgl['pore.coords'][pair[0]][0] for pair in connected_pores]
x2 = [sgl['pore.coords'][pair[1]][0] for pair in connected_pores]
same_x = [x - y == 0 for x, y in zip(x1,x2)]
factor = [s*.95 + (not s)*1 for s in same_x]
throat_diameters = sgl['throat.diameter'][throats]*factor
sgl.set_data(throats = throats, prop = 'diameter', data = throat_diameters)

#reset aspects relying on pore and throat sizes
geo.pore_volume()
geo.throat_volume()
geo.throat_surface_area()

#set up fluids 
air = OpenPNM.Fluids.Air(network = sgl, name = 'air')
water = OpenPNM.Fluids.Water(network = sgl, name = 'water')
#MYSTERIOUSLY BROKEN LINE
#water_sgl.add_property(prop = 'contact_angle', model = 'constant', value = 100)

sgl.regenerate_fluids()
water['pore.contact_angle'] = 100

#set up physics 
phys_water = OpenPNM.Physics.GenericPhysics(network=sgl,fluid=water, geometry = geo, name='standard_water_physics')
phys_air = OpenPNM.Physics.GenericPhysics(network=sgl,fluid=air, geometry = geo, name='standard_air_physics')
   
phys_water.add_property(prop='capillary_pressure', model='washburn')
phys_water.add_property(prop='hydraulic_conductance',model='hagen_poiseuille')
phys_water.add_property(prop='diffusive_conductance', model='bulk_diffusion', shape = 'square')

phys_air.add_property(prop='hydraulic_conductance',model='hagen_poiseuille')
phys_air.add_property(prop='diffusive_conductance', model='bulk_diffusion') 


sgl.regenerate_physics()

#late pore filling?!

#running invasion percolation
sat = sat = {'00': [], '10': [], '20': [], '01': [], '11': [], '21': []}
perm_air = {'00': [], '10': [], '20': [], '01': [], '11': [], '21': []}
diff_air = {'00': [], '10': [], '20': [], '01': [], '11': [], '21': []}
perm_water = {'00': [], '10': [], '20': [], '01': [], '11': [], '21': []}
diff_water = {'00': [], '10': [], '20': [], '01': [], '11': [], '21': []}

inlets = sgl.get_pore_indices(labels = ['bottom','boundary'],mode='intersection')
outlets = sgl.get_pore_indices(labels = ['top','boundary'],mode='intersection')


end_condition = 'breakthrough'
    
OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=sgl,loglevel=30)
OP_1.setup(invading_fluid = water, defending_fluid = air, inlets = inlets,npts=100)
OP_1.run() 


max_inv_seq = max(OP_1.get_throat_data(prop = 'inv_seq'))

for x in range(21):
    OP_1.update(seq = max_inv_seq*(x/20))
    
    print('seq = '+str(round(max_inv_seq*(x/20)))+' Pa out of '+str(round(max_inv_seq))+' total sequences')
    
    for mode_increment in range(len(modes)):
        #adding multiphase conductances
        phys_air.add_property(prop='multiphase',model='conduit_conductance',
                          conductance = 'diffusive_conductance', prop_name='conduit_diffusive_conductance',mode=modes[mode_increment])
        phys_water.add_property(prop='multiphase',model='conduit_conductance',
                          conductance = 'diffusive_conductance', prop_name='conduit_diffusive_conductance',mode=modes[mode_increment])
        phys_air.add_property(prop='multiphase',model='conduit_conductance',
                          conductance = 'hydraulic_conductance', prop_name='conduit_hydraulic_conductance',mode=modes[mode_increment])
        phys_water.add_property(prop='multiphase',model='conduit_conductance',
                          conductance = 'hydraulic_conductance', prop_name='conduit_hydraulic_conductance',mode=modes[mode_increment])
        sgl.regenerate_physics()
        
        
        #run Stokes Flow and find Permeability
        #single phase
        Stokes_alg_single_phase_air = OpenPNM.Algorithms.StokesFlow(loggername = 'Stokes', name = 'Stokes_alg_single_phase_air', network = sgl, loglevel = 40)
        Stokes_alg_single_phase_water = OpenPNM.Algorithms.StokesFlow(loggername = 'Stokes_2', name = 'Stokes_alg_single_phase_water', network = sgl, loglevel = 40)
        
        Fickian_alg_single_phase_air = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian', name = 'Fickian_alg_single_phase_air', network = sgl, loglevel = 40)
        Fickian_alg_single_phase_water = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian_2', name = 'Fickian_alg_single_phase_water', network = sgl, loglevel = 40)
        
        Stokes_alg_multi_phase_air = OpenPNM.Algorithms.StokesFlow(loggername = 'Stokes', name = 'Stokes_alg_multi_phase_air', network = sgl, loglevel = 40)
        Stokes_alg_multi_phase_water = OpenPNM.Algorithms.StokesFlow(loggername = 'Stokes_2', name = 'Stokes_alg_multi_phase_water', network = sgl, loglevel = 40)
        
        Fickian_alg_multi_phase_air = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian', name = 'Fickian_alg_multi_phase_air', network = sgl, loglevel = 40)
        Fickian_alg_multi_phase_water = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian_2', name = 'Fickian_alg_multi_phase_water', network = sgl, loglevel = 40)
        
        #setting up boundary conditions and calculating effective_permeability
        #BC1
        
        bounds = [['front', 'back'], ['left', 'right'], ['top', 'bottom']]
        for bound_increment in range(len(bounds)):
            BC1_pores = sgl.pores(labels=[bounds[bound_increment][0], 'boundary'],mode='intersection')
            BC2_pores = sgl.pores(labels=[bounds[bound_increment][1], 'boundary'],mode='intersection')
            
            
            Stokes_alg_single_phase_air.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
            Stokes_alg_single_phase_water.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
            Fickian_alg_single_phase_air.set_boundary_conditions(bctype = 'Dirichlet', bcvalue = .6, pores = BC1_pores)
            Fickian_alg_single_phase_water.set_boundary_conditions(bctype = 'Dirichlet', bcvalue = .6, pores = BC1_pores)
            
            Stokes_alg_multi_phase_air.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
            Stokes_alg_multi_phase_water.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
            Fickian_alg_multi_phase_air.set_boundary_conditions(bctype = 'Dirichlet', bcvalue = .6, pores = BC1_pores)
            Fickian_alg_multi_phase_water.set_boundary_conditions(bctype = 'Dirichlet', bcvalue = .6, pores = BC1_pores)   
            
            #BC2    
            Stokes_alg_single_phase_air.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.2, pores=BC2_pores)
            Stokes_alg_single_phase_water.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.2, pores=BC2_pores)
            Fickian_alg_single_phase_air.set_boundary_conditions(bctype = 'Dirichlet', bcvalue = .2, pores = BC2_pores)
            Fickian_alg_single_phase_water.set_boundary_conditions(bctype = 'Dirichlet', bcvalue = .2, pores = BC2_pores)
            
            Stokes_alg_multi_phase_air.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.2, pores=BC2_pores)
            Stokes_alg_multi_phase_water.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.2, pores=BC2_pores)
            Fickian_alg_multi_phase_air.set_boundary_conditions(bctype = 'Dirichlet', bcvalue = .2, pores = BC2_pores)
            Fickian_alg_multi_phase_water.set_boundary_conditions(bctype = 'Dirichlet', bcvalue = .2, pores = BC2_pores)
            
            #conduit conductance
            Stokes_alg_single_phase_air.setup(conductance = 'hydraulic_conductance',fluid=air)
            Stokes_alg_single_phase_water.setup(conductance = 'hydraulic_conductance',fluid=water)
            Fickian_alg_single_phase_air.setup(conductance = 'diffusive_conductance',fluid=air) 
            Fickian_alg_single_phase_water.setup(conductance = 'diffusive_conductance',fluid=water)
            
            Stokes_alg_multi_phase_air.setup(conductance = 'conduit_hydraulic_conductance',fluid=air)
            Stokes_alg_multi_phase_water.setup(conductance = 'conduit_hydraulic_conductance',fluid=water)
            Fickian_alg_multi_phase_air.setup(conductance = 'conduit_diffusive_conductance',fluid=air) 
            Fickian_alg_multi_phase_water.setup(conductance = 'conduit_diffusive_conductance',fluid=water)
            
            #run
            Stokes_alg_single_phase_air.run()
            Stokes_alg_single_phase_water.run()
            Fickian_alg_single_phase_air.run()
            Fickian_alg_single_phase_water.run()
            
            Stokes_alg_multi_phase_air.run()
            Stokes_alg_multi_phase_water.run()
            Fickian_alg_multi_phase_air.run()
            Fickian_alg_multi_phase_water.run()
            
            #calc effective properties
            effective_permeability_air_single = Stokes_alg_single_phase_air.calc_eff_permeability(clean = False)  
            effective_diffusivity_air_single = Fickian_alg_single_phase_air.calc_eff_diffusivity(clean = False)
            effective_permeability_water_single = Stokes_alg_single_phase_water.calc_eff_permeability(clean = False)  
            effective_diffusivity_water_single = Fickian_alg_single_phase_water.calc_eff_diffusivity(clean = False)
            
            effective_permeability_air_multi = Stokes_alg_multi_phase_air.calc_eff_permeability(clean = False)  
            effective_diffusivity_air_multi = Fickian_alg_multi_phase_air.calc_eff_diffusivity(clean = False)
            effective_permeability_water_multi = Stokes_alg_multi_phase_water.calc_eff_permeability(clean = False)  
            effective_diffusivity_water_multi = Fickian_alg_multi_phase_water.calc_eff_diffusivity(clean = False)
            
            final_pores = water.get_pore_data('occupancy')*1
            pore_volumes = sgl.get_pore_data(prop = 'volume')
            final_throats = water.get_throat_data('occupancy')*1
            throat_volumes = sgl.get_throat_data(prop = 'volume')
            
            saturation = (sum(final_pores*pore_volumes) + sum(final_throats*throat_volumes))/(sum(pore_volumes) + sum(throat_volumes))
            
            relative_eff_perm_air = effective_permeability_air_multi/effective_permeability_air_single
            relative_eff_perm_water = effective_permeability_water_multi/effective_permeability_water_single
            relative_eff_diff_air = effective_diffusivity_air_multi/effective_diffusivity_air_single
            relative_eff_diff_water = effective_diffusivity_water_multi/effective_diffusivity_water_single
            
            sat[str(bound_increment) + str(mode_increment)].append(saturation)
            perm_air[str(bound_increment) + str(mode_increment)].append(relative_eff_perm_air)
            diff_air[str(bound_increment) + str(mode_increment)].append(relative_eff_diff_air) 
            perm_water[str(bound_increment) + str(mode_increment)].append(relative_eff_perm_water)
            diff_water[str(bound_increment) + str(mode_increment)].append(relative_eff_diff_water)

from matplotlib.font_manager import FontProperties

fontP = FontProperties()
fontP.set_size('small')
#setting up subplots
#f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
fig = plt.figure(num=1, figsize=(10, 10), dpi=80, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(221)   #top left
ax2 = fig.add_subplot(222)   #top right
ax3 = fig.add_subplot(223)   #bottom left
ax4 = fig.add_subplot(224)

x_values1 = [x/20 for x in range(21)]
z = '.75'

#plots for subplot1 - loose permeability
p1, = ax1.plot(sat['00'], perm_water['00'], color = 'k', linestyle = '-', marker = 'o')
p2, = ax1.plot(sat['10'], perm_water['10'], color = z, linestyle = '-', marker = 'o')
p3, = ax1.plot(sat['20'], perm_water['20'], color = 'w', linestyle = '-', marker = 'o')
p4, = ax1.plot(sat['00'], perm_air['00'], color = 'k', linestyle = '-', marker = 'D')
p5, = ax1.plot(sat['10'], perm_air['10'], color = z, linestyle = '-', marker = 'D')
p6, = ax1.plot(sat['20'], perm_air['20'], color = 'w', linestyle = '-', marker = 'D')
p10, = ax1.plot(x_values1, [x**(3) for x in x_values1], 'k--')
ax1.plot(x_values1, [(1-x)**(3) for x in x_values1], 'k--')
ax1.set_title("loose permeability")
ax1.set_ylabel('permeability')
ax1.set_ylim([0,1])


#plots for subplot2 - strict permeability
p1, = ax2.plot(sat['00'], perm_water['00'], color = 'k', linestyle = '-', marker = 'o')
p2, = ax2.plot(sat['10'], perm_water['10'], color = z, linestyle = '-', marker = 'o')
p3, = ax2.plot(sat['20'], perm_water['20'], color = 'w', linestyle = '-', marker = 'o')
p7, = ax2.plot(sat['01'], perm_air['01'], color = 'k', linestyle = '-', marker = '^')
p8, = ax2.plot(sat['11'], perm_air['11'], color = z, linestyle = '-', marker = '^')
p9, = ax2.plot(sat['21'], perm_air['21'], color = 'w', linestyle = '-', marker = '^')
p10, = ax2.plot(x_values1, [x**(3) for x in x_values1], 'k--')
ax2.plot(x_values1, [(1-x)**(3) for x in x_values1], 'k--')
ax2.set_title("strict permeability")
ax2.set_ylim([0,1])

#need to work on legend to match up with the right things
lgd1 = ax2.legend([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10],
           ["KrWater,x", "KrWater,y", "KrWater,z",
           "KrAir,x (loose)","KrAir,y (loose)","KrAir,z (loose)",
           "KrAir,x (strict)","KrAir,y (strict)","KrAir,z (strict)", "a = 3"], loc='center left', bbox_to_anchor=(1, 0.5), prop = fontP)

#plots for subplot3 - loose diffusivity
p11, = ax3.plot(sat['00'], diff_water['00'], color = 'k', linestyle = '-', marker = 'o')
p12, = ax3.plot(sat['10'], diff_water['10'], color = z, linestyle = '-', marker = 'o')
p13, = ax3.plot(sat['20'], diff_water['20'], color = 'w', linestyle = '-', marker = 'o')
p14, = ax3.plot(sat['00'], diff_air['00'], color = 'k', linestyle = '-', marker = 'D')
p15, = ax3.plot(sat['10'], diff_air['10'], color = z, linestyle = '-', marker = 'D')
p16, = ax3.plot(sat['20'], diff_air['20'], color = 'w', linestyle = '-', marker = 'D')
p20, = ax3.plot(x_values1, [x**(2) for x in x_values1], 'k--')
ax3.plot(x_values1, [(1-x)**(2) for x in x_values1], 'k--')
ax3.set_title("loose diffusivity")
ax3.set_ylabel("diffusivity")
ax3.set_xlabel("saturation")
ax3.set_ylim([0,1])


#plots for subplot4 - strict diffusivity
p11, = ax4.plot(sat['00'], diff_water['00'], color = 'k', linestyle = '-', marker = 'o')
p12, = ax4.plot(sat['10'], diff_water['10'], color = z, linestyle = '-', marker = 'o')
p13, = ax4.plot(sat['20'], diff_water['20'], color = 'w', linestyle = '-', marker = 'o')
p17, = ax4.plot(sat['01'], diff_air['01'], color = 'k', linestyle = '-', marker = '^')
p18, = ax4.plot(sat['11'], diff_air['11'], color = z, linestyle = '-', marker = '^')
p19, = ax4.plot(sat['21'], diff_air['21'], color = 'w', linestyle = '-', marker = '^')
p20, = ax4.plot(x_values1, [x**(2) for x in x_values1], 'k--')
ax4.plot(x_values1, [(1-x)**(2) for x in x_values1], 'k--')
ax4.set_title("strict diffusivity")
ax4.set_xlabel("saturation")
ax4.set_ylim([0,1])

lgd2 = ax4.legend([p11, p12, p13, p14, p15, p16, p17, p18, p19, p20],
           ["DrWater,x", "DrWater,y", "DrWater,z",
           "DrAir,x (loose)","DrAir,y (loose)","DrAir,z (loose)",
           "DrAir,x (strict)","DrAir,y (strict)","DrAir,z (strict)", "a = 3"], loc='center left', bbox_to_anchor=(1, 0.5), prop = fontP)

fig.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1) 
         
fig.show()


