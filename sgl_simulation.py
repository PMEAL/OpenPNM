# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 16:38:29 2014

@author: Jackie
"""
import OpenPNM
import matplotlib.pyplot as plt
import scipy as sp

def run_simulation():
    r"""
    returns array
    array[0] = saturation
    array[1] = relative permeability
    array[2] = relative effective diffusivity
    """  
    Lc = 40.5e-6
    N = 10
    
    #setting up network
    sgl = OpenPNM.Network.Cubic(name = 'SGL10BA')
    sgl.generate(divisions = [N, N, N], add_boundaries = False, lattice_spacing = [Lc])
    
    #pore size distribution parameters
    lmbda = 9e-6
    k = 3.5
    bmin = 9e-6
    xmax = .9
    
    #set up geometries
    geo = OpenPNM.Geometry.GenericGeometry(name = 'sgl', network = sgl)
    geo.set_locations(pores = sgl.pores('all'), throats = 'all')
    
    Np = sgl.num_pores('sgl')
    value = lmbda*(-sp.log(1-sp.random.rand(Np)*xmax))**(-1/k) + bmin
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
    #phys_water.add_property(prop='multiphase',model='conduit_conductance',
                      #conductance = 'diffusive_conductance', prop_name='conduit_diffusive_conductance',mode='strict')
    
    phys_air.add_property(prop='hydraulic_conductance',model='hagen_poiseuille')
    phys_air.add_property(prop='diffusive_conductance', model='bulk_diffusion') 
    #phys_air.add_property(prop='multiphase',model='conduit_conductance',
                      #conductance = 'diffusive_conductance', prop_name='conduit_diffusive_conductance',mode='strict')

    
    sgl.regenerate_physics()
    
    #late pore filling?!
    
    #running invasion percolation
    inlets = sgl.get_pore_indices(labels = ['bottom'])
    outlets = sgl.get_pore_indices(labels = ['top'])
    
    end_condition = 'breakthrough'
    
    OP_1 = OpenPNM.Algorithms.InvasionPercolation(network = sgl, name = 'OP_1',loglevel=30)
    OP_1.setup(invading_fluid = water, defending_fluid = air, inlets = inlets, outlets = outlets, end_condition = end_condition)
    OP_1.run()    
    OP_1.update()
    
    #run Stokes Flow and find Permeability
    Stokes_alg = OpenPNM.Algorithms.StokesFlow(loggername = 'Stokes', name = 'Stokes_alg', network = sgl)
    Stokes_alg_2 = OpenPNM.Algorithms.StokesFlow(loggername = 'Stokes_2', name = 'Stokes_alg_2', network = sgl)
    
    Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian', name = 'fickian_alg', network = sgl)
    Fickian_alg_2 = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian_2', name = 'fickian_alg_2', network = sgl)
    
    #setting up boundary conditions and calculating effective_permeability
    BC1_pores = sgl.pores(labels=['top','front'],mode='intersection')
    Stokes_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
    Stokes_alg_2.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
    Fickian_alg.set_boundary_conditions(bctype = 'Dirichlet', bcvalue = .6, pores = BC1_pores)
    Fickian_alg_2.set_boundary_conditions(bctype = 'Dirichlet', bcvalue = .6, pores = BC1_pores)   
    
    BC2_pores = sgl.pores(labels=['top','back'],mode='intersection')
    Stokes_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.2, pores=BC2_pores)
    Stokes_alg_2.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.2, pores=BC2_pores)
    Fickian_alg.set_boundary_conditions(bctype = 'Dirichlet', bcvalue = .2, pores = BC2_pores)
    Fickian_alg_2.set_boundary_conditions(bctype = 'Dirichlet', bcvalue = .2, pores = BC2_pores)
    
    Stokes_alg.setup(fluid = air)
    Stokes_alg_2.setup(fluid = water)
    Fickian_alg.setup(fluid = air) 
    Fickian_alg_2.setup(fluid = water)
    Stokes_alg.run()
    Stokes_alg_2.run()
    Fickian_alg.run()
    Fickian_alg_2.run()
    
    effective_permeability = Stokes_alg.calc_eff_permeability(clean = False)  
    effective_diffusivity = Fickian_alg.calc_eff_diffusivity(clean = False)
    effective_permeability_water = Stokes_alg_2.calc_eff_permeability(clean = False)  
    effective_diffusivity_water = Fickian_alg_2.calc_eff_diffusivity(clean = False)
    
    final_pores = water.get_pore_data('occupancy')*1
    pore_volumes = sgl.get_pore_data(prop = 'volume')
    
    sum_volume = 0
    filled_volume = 0
    
    for i in range(len(sgl.pores())):
        sum_volume += pore_volumes[i]
        if final_pores[i] != 0:
            filled_volume += pore_volumes[i]
        
    saturation = filled_volume/sum_volume
    
    OP_1.update(IPseq = 0)
    Stokes_alg.setup(fluid = air)
    Stokes_alg_2.setup(fluid = water)
    Fickian_alg.setup(fluid = air) 
    Fickian_alg_2.setup(fluid = water)
    Stokes_alg.run()
    Fickian_alg.run()
    Stokes_alg_2.run()
    Fickian_alg_2.run()
    
    effective_permeability_2 = Stokes_alg.calc_eff_permeability(clean = False)
    effective_diffusivity_2 = Fickian_alg.calc_eff_diffusivity(clean = False)
    effective_permeability_water_2 = Stokes_alg_2.calc_eff_permeability(clean = False)  
    effective_diffusivity_water_2 = Fickian_alg_2.calc_eff_diffusivity(clean = False)
    
    return[saturation, effective_permeability/effective_permeability_2[0], effective_permeability_water/effective_permeability_water_2[0], effective_diffusivity/effective_diffusivity_2, effective_diffusivity_water/effective_diffusivity_water_2]
    
def saving_simulation_data(number_times_run, x_values, y_values, y2_values, z_values, z2_values):
    
    for x in range(number_times_run):
        data = run_simulation()
        x_values.append(data[0])
        y_values.append(data[1][0])
        y2_values.append(data[2][0])
        z_values.append(data[3][0])
        z2_values.append(data[4][0])
        
        if(x == 0):
            f = open('simulationResults', 'w')
        else:
            f = open('simulationResults', 'a')
        f.write('\nsimulation: ')
        f.write(str(x))
        f.write('\nsaturation: ')
        f.write(str(data[0]))
        f.write('\relative permeability: ')
        f.write(str(data[1]))
        f.close()
        
x_values = []
y_values = []
y2_values = []
z_values = []
z2_values = []
saving_simulation_data(number_times_run = 1, x_values = x_values, y_values = y_values, y2_values = y2_values, z_values = z_values, z2_values = z2_values)

x_values1 = [x/20 for x in range(20)]

plt.plot(x_values, y_values, 'ro')
plt.plot(x_values, y2_values, 'bs')
plt.plot(x_values1, [x**(3) for x in x_values1], 'g--')

plt.title('saturation versus relative permeability')
plt.xlabel('saturation')
plt.ylabel('relative permeability')
plt.show()

plt.plot(x_values, z_values, 'ro')
plt.plot(x_values, z2_values, 'bs')
plt.plot(x_values1, [x**(2) for x in x_values1], 'g--')

plt.title('saturation versus relative diffusivity')
plt.xlabel('saturation')
plt.ylabel('relative diffusivity')
plt.show()