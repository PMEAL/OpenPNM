# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 16:38:29 2014

@author: Jackie
"""
import OpenPNM
import matplotlib.pyplot as plt
import scipy as sp


r"""
returns array
array[0] = saturation
array[1] = relative permeability
array[2] = relative effective diffusivity
"""  
Lc = 40.5e-6
N = 10


#setting up network
sgl = OpenPNM.Network.Cubic(name = 'SGL10BA', loglevel = 40)
sgl.generate(divisions = [N, N, N], add_boundaries = True, lattice_spacing = [Lc])

#pore size distribution parameters
lmbda = 9e-6
k = 3.5
bmin = 9e-6
xmax = .9

#set up geometries
geo = OpenPNM.Geometry.GenericGeometry(name = 'sgl', network = sgl)
geo.set_locations(pores = sgl.pores('all'), throats = 'all')

Np = sgl.num_pores('sgl')
value = [min(lmbda*(-sp.log(1-sp.random.rand()*xmax))**(-1/k) + bmin, Lc) for x in sgl.pores()]
sgl.set_data(prop='diameter',pores=geo.pores(),data=value)

geo.add_method(prop='pore_volume',model='sphere')
geo.add_method(prop='throat_length',model='straight')
geo.add_method(prop='throat_volume',model='cuboid')

connected_pores = sgl.find_connected_pores(geo.throats())
value = [min(sgl.get_data(pores = pair[0], prop = 'diameter'), sgl.get_data(pores = pair[1], prop = 'diameter')) for pair in connected_pores]
sgl.set_data(prop = 'diameter', throats = geo.throats(), data = value)


sgl.regenerate_geometries()

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

#set up fluids 
air = OpenPNM.Fluids.Air(network = sgl, name = 'air')
water = OpenPNM.Fluids.Water(network = sgl, name = 'water')
#MYSTERIOUSLY BROKEN LINE
#water_sgl.add_property(prop = 'contact_angle', model = 'constant', value = 100)
water['pore.contact_angle'] = 100
sgl.regenerate_fluids()

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


max_capillary_pressure = max(OP_1.get_throat_data(prop = 'inv_Pc'))
#    OP_1.update(max(OP_1.get_throat_data(prop='inv_Pc')))
#    OP_1.update(sp.stats.mstats.mode(OP_1.get_pore_data(prop='inv_Pc')))
    
for x in range(100):
    OP_1.update(max_capillary_pressure*(x/100.0))
    print('Pc = '+str(round(max_capillary_pressure*(x/100.0)))+' Pa out of '+str(round(max_capillary_pressure))+' Pa')
    modes = ['loose', 'strict']
    
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
            
            sum_volume = 0
            filled_volume = 0
            
            for i in range(len(sgl.pores())):
                sum_volume += pore_volumes[i]
                if final_pores[i] != 0:
                    filled_volume += pore_volumes[i]
                    
            for j in range(len(sgl.throats())):
                sum_volume += throat_volumes[j]
                if final_throats[j] != 0:
                    filled_volume += throat_volumes[j]
                
            saturation = filled_volume/sum_volume
            
            relative_eff_perm_air = effective_permeability_air_multi/effective_permeability_air_single
            relative_eff_perm_water = effective_permeability_water_multi/effective_permeability_water_single
            relative_eff_diff_air = effective_diffusivity_air_multi/effective_diffusivity_air_single
            relative_eff_diff_water = effective_diffusivity_water_multi/effective_diffusivity_water_single
            
            sat[str(bound_increment) + str(mode_increment)].append(saturation)
            perm_air[str(bound_increment) + str(mode_increment)].append(relative_eff_perm_air)
            diff_air[str(bound_increment) + str(mode_increment)].append(relative_eff_diff_air) 
            perm_water[str(bound_increment) + str(mode_increment)].append(relative_eff_perm_water)
            diff_water[str(bound_increment) + str(mode_increment)].append(relative_eff_diff_water)

x_values1 = [x/20 for x in range(21)]
z = '.75'

p1, = plt.plot(sat['00'], perm_water['00'], color = 'k', linestyle = '-', marker = 'o')
p2, = plt.plot(sat['10'], perm_water['10'], color = z, linestyle = '-', marker = 'o')
p3, = plt.plot(sat['20'], perm_water['20'], color = 'w', linestyle = '-', marker = 'o')
p4, = plt.plot(sat['00'], perm_air['00'], color = 'k', linestyle = '-', marker = 'D')
p5, = plt.plot(sat['10'], perm_air['10'], color = z, linestyle = '-', marker = 'D')
p6, = plt.plot(sat['20'], perm_air['20'], color = 'w', linestyle = '-', marker = 'D')
p7, = plt.plot(sat['01'], perm_air['01'], color = 'k', linestyle = '-', marker = '^')
p8, = plt.plot(sat['11'], perm_air['11'], color = z, linestyle = '-', marker = '^')
p9, = plt.plot(sat['21'], perm_air['21'], color = 'w', linestyle = '-', marker = '^')
p10, = plt.plot(x_values1, [x**(3) for x in x_values1], 'k--')
plt.plot(x_values1, [(1-x)**(3) for x in x_values1], 'k--')

plt.legend([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10],
           ["KrWater,x", "KrWater,y", "KrWater,z",
           "KrAir,x (loose)","KrAir,y (loose)","KrAir,z (loose)",
           "KrAir,x (strict)","KrAir,y (strict)","KrAir,z (strict)", "a = 3"], 
loc = 'lower left', bbox_to_anchor = (1, .05)) 


#plt.title('relative permeability versus saturation')
plt.xlabel('saturation')
plt.ylabel('relative permeability')
plt.show()

p1, = plt.plot(sat['00'], diff_water['00'], color = 'k', linestyle = '-', marker = 'o')
p2, = plt.plot(sat['10'], diff_water['10'], color = z, linestyle = '-', marker = 'o')
p3, = plt.plot(sat['20'], diff_water['20'], color = 'w', linestyle = '-', marker = 'o')
p4, = plt.plot(sat['00'], diff_air['00'], color = 'k', linestyle = '-', marker = 'D')
p5, = plt.plot(sat['10'], diff_air['10'], color = z, linestyle = '-', marker = 'D')
p6, = plt.plot(sat['20'], diff_air['20'], color = 'w', linestyle = '-', marker = 'D')
p7, = plt.plot(sat['01'], diff_air['01'], color = 'k', linestyle = '-', marker = '^')
p8, = plt.plot(sat['11'], diff_air['11'], color = z, linestyle = '-', marker = '^')
p9, = plt.plot(sat['21'], diff_air['21'], color = 'w', linestyle = '-', marker = '^')
p10, = plt.plot(x_values1, [x**(2) for x in x_values1], 'k--')
plt.plot(x_values1, [(1-x)**(2) for x in x_values1], 'k--')

plt.legend([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10],
           ["DrWater,x", "DrWater,y", "DrWater,z",
           "DrAir,x (loose)","DrAir,y (loose)","DrAir,z (loose)",
           "DrAir,x (strict)","DrAir,y (strict)","DrAir,z (strict)", "a = 2"],
           loc = 'lower left', bbox_to_anchor = (1, .05))

#plt.title('relative diffusivity versus saturation')
plt.xlabel('saturation')
plt.ylabel('relative diffusivity')
plt.show()