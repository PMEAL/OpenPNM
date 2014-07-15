
# coding: utf-8

## Regerating Data from Gostick et Al, "Pore network modeling of fibrous gas diffusion layers for polymer electrolyte membrane fule cells"

# In this tutorial, we will regenerate data from Gostick's 2007 paper.  This will both show that OpenPNM can recreate results accurately, and will also show some more specific uses of OpenPNM.  While his paper deals with both SGL and Toray GDLs, we will deal only with SGL. 
# 
# There will be a general layout to complete this simulation:
# 1) Set up network
# 2) Set up geometry and geometrical methods
# 3) constrict throat's by a constriction factor
# 4) Set up fluids and methods
# 5) Set up fluid physics and methods
# 6) Run invasion percolation
# 7) Run Stokes and Fickian algorithms 
# 8) generate effective permeability and effective diffusivity values at different saturations
# 9) plot generated data
# 
# To begin, we must download OpenPNM found at https://github.com/PMEAL/OpenPNM.  Then we can begin our simulation by adding OpenPNM to our pythonpath, and importing OpenPNM.  
# 
# The command line argument for adding OpenPNM to our python path (for MacOSX, note that a folder along the path used must contain an init folder in order for this to work):
# export PYTHONPATH=$PYTHONPATH:/path/to/OpenPNM
# 
# We should also import matplotlib.pyplot so that we can plot our graphs at the end.

# In[1]:

cd OpenPNM


# In[2]:

import OpenPNM
import matplotlib.pyplot as plt


### Setting up Network and Geometry

# To begin our simulation, we must first generate our SGL network and geometry.  This includes 
# 
# 1) creating a cubic network object and an SGL10 geometry object
# 2) sending our geometry object our internal pores
# 3) calculating values for throat and pore properties for both internal and boundary pores
# 4) accounting for pores that are too big (making maximum pore size the lattice parameter)

# In[3]:

Lc = 40.5e-6


#setting up network
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


# Before we move on to setting up our fluid and physics objects, we must constrict throats in the z and y direction.  For his SGL simulation, he uses a constriction factor of .95.  Finally, because we have changed values for pore and throat diameters (first by accounting for pores and throats that are too big, and the finally constricting throats in the y and z directions), we must recalculate all pore and throat values relying on these diameters.

# In[4]:

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
geo.throat_length()
geo.pore_volume()
geo.throat_volume()
geo.throat_surface_area()


### Setting up fluids and physics

# Now we are ready to set up our fluids (water and air) and the physics corresponding to each of these fluids.  OpenPNM has built in air and water objects, so we can use those.  We must remember to use regenerate_fluids() so that important values for each fluid will be set.  However, Gostick specifies using a water pore contact angle of 100, so we will reset this value after regenerating our fluids.

# In[5]:

#set up fluids 
air = OpenPNM.Fluids.Air(network = sgl, name = 'air')
water = OpenPNM.Fluids.Water(network = sgl, name = 'water')

sgl.regenerate_fluids()
water['pore.contact_angle'] = 100


# We are now ready to establish physical properties for our fluid objects.  To do this, we will:
# 1) create physics objects associated with our fluids
# 2) add methods concerning each physical property we will want calculated
# 3) use our regenerate_physics() method to recalculate these properties

# In[6]:

#set up physics 
phys_water = OpenPNM.Physics.GenericPhysics(network=sgl,fluid=water, geometry = geo, name='standard_water_physics')
phys_air = OpenPNM.Physics.GenericPhysics(network=sgl,fluid=air, geometry = geo, name='standard_air_physics')
   
phys_water.add_property(prop='capillary_pressure', model='washburn')
phys_water.add_property(prop='hydraulic_conductance',model='hagen_poiseuille')
phys_water.add_property(prop='diffusive_conductance', model='bulk_diffusion', shape = 'square')

phys_air.add_property(prop='hydraulic_conductance',model='hagen_poiseuille')
phys_air.add_property(prop='diffusive_conductance', model='bulk_diffusion') 


sgl.regenerate_physics()


### Running ordinary percolation, Fickian diffusion, and Stokes flow

# Gostick uses ordinary percolation to spread water through his GDL before calculating relative permeability and relative diffusivity.  This way, a graph showing the relationship between saturation and relative permeability and between saturation and relative diffusivity can be created.  
# 
# To run our ordinary percolation, we will:
# 1) pick inlet and outlet pores
# 2) create an Ordinary Percolation algorithm object
# 3) setup our algorithm object
# 4) run our algorithm object
# 5) call update() so that occupancy of pores and throats for each fluid will be set

# In[7]:

inlets = sgl.get_pore_indices(labels = ['bottom','boundary'],mode='intersection')

#using every other pore in the bottom and boundary as an inlet to prevent extremely small values in the z direction
used_inlets = [inlets[x] for x in range(0, len(inlets), 2)]
    
OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=sgl,loglevel=30)
OP_1.setup(invading_fluid = water, defending_fluid = air, inlets = used_inlets,npts=100)
OP_1.run() 
OP_1.update()


# update() takes either a sequence value or a pressure value.  If we give no parameters, it uses a capillary pressure value of 0.  This means that the percolation process is updated to the point in time at which the capillary pressure is 0, or when all pores have not been invaded by the invading fluid.  We can check that water pore occupancy is 0 for all pores.

# In[8]:

water['pore.occupancy'][water.pores()]


# The next step will be to calculate effective diffusivity and permeability at different saturations.  Note that we want to run Fickian diffusion and Stokes flow algorithms at different points within our ordinary percolation process.  OpenPNM has a very helpful update() method for updating the occupancy of pores to their values during a specified part of percolation.  During percolation, each pore is given a sequence value showing when in time it was invaded.  We can send update() a sequence parameter, determining when during the percolation we want to update our pore occupancy to.  
# 
# The rest of our code will exist within a loop updating our network to different stages of percolation, so that we may view our relative diffusivity and permeability at different points of saturation.
# 
# Before we add in the loop aspect, we will walk through the code that will be inside the loop.  
# 
# First, we will want to add a physics property that recalculates diffusive and hydraulic conductance in each throat based on occupancy.  There are two modes for this- 'strict' and 'loose'.  'loose' lowers conductivity if the throat is filled with invading fluid, while 'strict' lowers conductivity if the throat or either neighboring pore is filled with the invading fluid.  We will start off with 'loose', and later try 'strict' to compare results.

# In[9]:

#adding multiphase conductances
phys_air.add_property(prop='multiphase',model='conduit_conductance',
                  conductance = 'diffusive_conductance', prop_name='conduit_diffusive_conductance',mode='loose')
phys_water.add_property(prop='multiphase',model='conduit_conductance',
                  conductance = 'diffusive_conductance', prop_name='conduit_diffusive_conductance',mode='loose')
phys_air.add_property(prop='multiphase',model='conduit_conductance',
                  conductance = 'hydraulic_conductance', prop_name='conduit_hydraulic_conductance',mode='loose')
phys_water.add_property(prop='multiphase',model='conduit_conductance',
                  conductance = 'hydraulic_conductance', prop_name='conduit_hydraulic_conductance',mode='loose')
sgl.regenerate_physics()


# We can finally instatiate, setup, and run our algorithm objects for Stokes flow and Fickian diffusion.  We want to set up 8 different algorithm objects.
# 
# 1) Stokes flow, single phase air
# 2) Stokes flow, multi phase air 
# 3) Stokes flow, single phase water
# 4) Stokes flow, multi phase water
# 5) Fickian diffusion, single phase air
# 6) Fickian diffusion, multi phase air 
# 7) Fickian diffusion, sing phase water
# 8) Fickian diffusion, multi phase water
# 
# Note that we want the algorithms that are single phase (where only the specified fluid exists in the network) to make our permeability and diffusivity values relative.  Any algorithm that is single phase will use the hydraulic or diffusive conductances before we recalculated based on occupancy.  This calls for our conductance parameter to be 'hydraulic_conductance' or 'diffusive_conductance' instead of 'conduit_hydraulic_conductance' or 'conduit_diffusive_conductance'.

# In[10]:

#run Stokes Flow and find Permeability
Stokes_alg_single_phase_air = OpenPNM.Algorithms.StokesFlow(loggername = 'Stokes', name = 'Stokes_alg_single_phase_air', network = sgl, loglevel = 30)
Stokes_alg_single_phase_water = OpenPNM.Algorithms.StokesFlow(loggername = 'Stokes_2', name = 'Stokes_alg_single_phase_water', network = sgl, loglevel = 30)

Fickian_alg_single_phase_air = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian', name = 'Fickian_alg_single_phase_air', network = sgl, loglevel = 30)
Fickian_alg_single_phase_water = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian_2', name = 'Fickian_alg_single_phase_water', network = sgl, loglevel = 30)

Stokes_alg_multi_phase_air = OpenPNM.Algorithms.StokesFlow(loggername = 'Stokes', name = 'Stokes_alg_multi_phase_air', network = sgl, loglevel = 30)
Stokes_alg_multi_phase_water = OpenPNM.Algorithms.StokesFlow(loggername = 'Stokes_2', name = 'Stokes_alg_multi_phase_water', network = sgl, loglevel = 30)

Fickian_alg_multi_phase_air = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian', name = 'Fickian_alg_multi_phase_air', network = sgl, loglevel = 30)
Fickian_alg_multi_phase_water = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian_2', name = 'Fickian_alg_multi_phase_water', network = sgl, loglevel = 30)

#setting up boundary conditions (diffusion will take place only in the Z direction for now)
BC1_pores = sgl.pores(labels=['top', 'boundary'],mode='intersection')
BC2_pores = sgl.pores(labels=['bottom', 'boundary'],mode='intersection')

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
Stokes_alg_single_phase_air.run(loglevel = 30)
Stokes_alg_single_phase_water.run(loglevel = 30)
Fickian_alg_single_phase_air.run(loglevel = 30)
Fickian_alg_single_phase_water.run(loglevel = 30)

Stokes_alg_multi_phase_air.run(loglevel = 30)
Stokes_alg_multi_phase_water.run(loglevel = 30)
Fickian_alg_multi_phase_air.run(loglevel = 30)
Fickian_alg_multi_phase_water.run(loglevel = 30)


# Now that we have run our all algorithms, we are ready to calculate saturation, all our effective properties, and finally our relative effective properties.  Because we called update() without setting any parameters, our GDL is filled with air.  Therefore, as a sanity check, any relative effective properties for air should be 1, while those for water should be 0.

# In[11]:

#saturation calculation
final_pores = water.get_pore_data('occupancy')*1
pore_volumes = sgl.get_pore_data(prop = 'volume')
final_throats = water.get_throat_data('occupancy')*1
throat_volumes = sgl.get_throat_data(prop = 'volume')

saturation = (sum(final_pores*pore_volumes) + sum(final_throats*throat_volumes))/(sum(pore_volumes) + sum(throat_volumes))

#calc effective properties
effective_permeability_air_single = Stokes_alg_single_phase_air.calc_eff_permeability(clean = False)  
effective_diffusivity_air_single = Fickian_alg_single_phase_air.calc_eff_diffusivity(clean = False)
effective_permeability_water_single = Stokes_alg_single_phase_water.calc_eff_permeability(clean = False)  
effective_diffusivity_water_single = Fickian_alg_single_phase_water.calc_eff_diffusivity(clean = False)

effective_permeability_air_multi = Stokes_alg_multi_phase_air.calc_eff_permeability(clean = False)  
effective_diffusivity_air_multi = Fickian_alg_multi_phase_air.calc_eff_diffusivity(clean = False)
effective_permeability_water_multi = Stokes_alg_multi_phase_water.calc_eff_permeability(clean = False)  
effective_diffusivity_water_multi = Fickian_alg_multi_phase_water.calc_eff_diffusivity(clean = False)

#calc relative effective properties
relative_eff_perm_air = effective_permeability_air_multi/effective_permeability_air_single
relative_eff_perm_water = effective_permeability_water_multi/effective_permeability_water_single
relative_eff_diff_air = effective_diffusivity_air_multi/effective_diffusivity_air_single
relative_eff_diff_water = effective_diffusivity_water_multi/effective_diffusivity_water_single

#printing values as a sanity check
print(relative_eff_perm_air)
print(relative_eff_perm_water)
print(relative_eff_diff_air)
print(relative_eff_diff_water)


# Now we will use the for loop above mentioned to find relative effective properties at different levels of saturation.  This way, we can update the network to different states of occupancy, save our saturation and relative effective property values into arrays, so that we can finally spit out a graph showing the trend.  To simplify and make our code run faster, we will show only 10 points along each curve, and only calculate effective properties in the Z direction.
# 
# Note that if we use the network size that Gostick has chosen (26x26x10) this may take a couple minutes.

# In[12]:

#initiating arrays for saving data gathered within loop
sat = []
perm_air = []
diff_air = []
perm_water = []
diff_water = []

max_inv_seq = max(OP_1.get_throat_data(prop = 'inv_seq'))

for x in range(11):
    OP_1.update(seq = max_inv_seq*(x/10))
    
    print('seq = '+str(round(max_inv_seq*(x/10)))+' Seq out of '+str(round(max_inv_seq))+' total sequences')

    #adding multiphase conductances
    phys_air.add_property(prop='multiphase',model='conduit_conductance',
                      conductance = 'diffusive_conductance', prop_name='conduit_diffusive_conductance',mode='loose')
    phys_water.add_property(prop='multiphase',model='conduit_conductance',
                      conductance = 'diffusive_conductance', prop_name='conduit_diffusive_conductance',mode='loose')
    phys_air.add_property(prop='multiphase',model='conduit_conductance',
                      conductance = 'hydraulic_conductance', prop_name='conduit_hydraulic_conductance',mode='loose')
    phys_water.add_property(prop='multiphase',model='conduit_conductance',
                      conductance = 'hydraulic_conductance', prop_name='conduit_hydraulic_conductance',mode='loose')
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
    
    BC1_pores = sgl.pores(labels=['top', 'boundary'],mode='intersection')
    BC2_pores = sgl.pores(labels=['bottom', 'boundary'],mode='intersection')


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

    sat.append(saturation)
    perm_air.append(relative_eff_perm_air)
    diff_air.append(relative_eff_diff_air) 
    perm_water.append(relative_eff_perm_water)
    diff_water.append(relative_eff_diff_water)


### Displaying final results

# We have all our data stored in our arrays, and we are ready to generate a graph.  This graph will be slightly less involved than Gosticks, as it will only include relative values in the Z directions, and only at 10 different points during the percolation.  For comparison purposes, we will setup our graphs to be 2 subplots showing each relative property for air and water.  On the figure will be a curve for air, a curve for water, and the curve that shows a general trend (saturation^2 or saturation^3).  A legend will help guide the eye.

# In[13]:

from matplotlib.font_manager import FontProperties

#setting font size small so that legends will fit property
fontP = FontProperties()
fontP.set_size('small')

#setting up subplots
#f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
fig = plt.figure(num=1, figsize=(10, 10), dpi=80, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(221)   #top left
ax2 = fig.add_subplot(223)   #bottom left

x_values1 = [x/20 for x in range(21)]
z = '.75'

#plots for subplot1 - loose permeability
p1, = ax1.plot(sat, perm_water, color = 'k', linestyle = '-', marker = 'o')
p2, = ax1.plot(sat, perm_air, color = 'k', linestyle = '-', marker = 'D')
p3, = ax1.plot(x_values1, [x**(3) for x in x_values1], 'k--')
ax1.plot(x_values1, [(1-x)**(3) for x in x_values1], 'k--')
ax1.set_title("loose permeability")
ax1.set_ylabel('permeability')
ax1.set_ylim([0,1])
ax1.set_xlim([0,1])

#need to work on legend to match up with the right things
lgd1 = ax1.legend([p1, p2, p3],
           ["KrWater", "KrAir (loose)", "a = 3"], loc='center left', bbox_to_anchor=(1, 0.5), prop = fontP)

#plots for subplot3 - loose diffusivity
p4, = ax2.plot(sat, diff_water, color = 'k', linestyle = '-', marker = 'o')
p5, = ax2.plot(sat, diff_air, color = 'k', linestyle = '-', marker = 'D')
p6, = ax2.plot(x_values1, [x**(2) for x in x_values1], 'k--')
ax2.plot(x_values1, [(1-x)**(2) for x in x_values1], 'k--')
ax2.set_title("loose diffusivity")
ax2.set_ylabel("diffusivity")
ax2.set_xlabel("saturation")
ax2.set_ylim([0,1])
ax2.set_xlim([0,1])

lgd2 = ax2.legend([p4, p5, p6],
           ["DrWater", "DrAir, (loose)", "a = 3"], loc='center left', bbox_to_anchor=(1, 0.5), prop = fontP)

fig.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1) 
         
fig.show()


# The figures generated thus far should appear something like this:

# In[2]:

from IPython.display import display
from IPython.display import Image

i = Image(url = 'http://i.imgur.com/wXz7xg9.png')
display(i)
g = Image(url = 'http://i.imgur.com/Gc76cz6.png')
display(g)


### 'Loose' and 'Strict' modes

# As mentioned above, OpenPNM has 'Loose' and 'Strict' modes for updating hydraulic and diffusive conductance based on occupancy.  To reiterate, 'Loose' will update the conductance when the throat is filled with water.  'Strict' will update the conductance when either the throat is filled with water, or at least one of the neighboring pores is filled with water.  We will run our for loop again, to gather data using 'Strict'.  A similar graph to before will be generated, for comparison purposes.

# In[14]:

#initiating arrays for saving data gathered within loop
sat = []
perm_air = []
diff_air = []
perm_water = []
diff_water = []

max_inv_seq = max(OP_1.get_throat_data(prop = 'inv_seq'))

for x in range(11):
    OP_1.update(seq = max_inv_seq*(x/10))
    
    print('seq = '+str(round(max_inv_seq*(x/10)))+' Seq out of '+str(round(max_inv_seq))+' total sequences')

    #adding multiphase conductances
    phys_air.add_property(prop='multiphase',model='conduit_conductance',
                      conductance = 'diffusive_conductance', prop_name='conduit_diffusive_conductance',mode='strict')
    phys_water.add_property(prop='multiphase',model='conduit_conductance',
                      conductance = 'diffusive_conductance', prop_name='conduit_diffusive_conductance',mode='strict')
    phys_air.add_property(prop='multiphase',model='conduit_conductance',
                      conductance = 'hydraulic_conductance', prop_name='conduit_hydraulic_conductance',mode='strict')
    phys_water.add_property(prop='multiphase',model='conduit_conductance',
                      conductance = 'hydraulic_conductance', prop_name='conduit_hydraulic_conductance',mode='strict')
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
    
    BC1_pores = sgl.pores(labels=['top', 'boundary'],mode='intersection')
    BC2_pores = sgl.pores(labels=['bottom', 'boundary'],mode='intersection')


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

    sat.append(saturation)
    perm_air.append(relative_eff_perm_air)
    diff_air.append(relative_eff_diff_air) 
    perm_water.append(relative_eff_perm_water)
    diff_water.append(relative_eff_diff_water)


# Finally, we will take this code generated using 'strict' mode and create identical graphs to before.  Note that the effective properties have a lower curve to them.  For example, relative effective permeability and diffusivity for air reach 0 at a lower saturation.  This makes sense, as it is harder for air to find a path if conductance is lowered even if one of the neighboring pores is filled with water.

# In[16]:

from matplotlib.font_manager import FontProperties

#setting font size small so that legends will fit property
fontP = FontProperties()
fontP.set_size('small')

#setting up subplots
#f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
fig = plt.figure(num=1, figsize=(10, 10), dpi=80, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(221)   #top left
ax2 = fig.add_subplot(223)   #bottom left

x_values1 = [x/20 for x in range(21)]
z = '.75'

#plots for subplot1 - loose permeability
p1, = ax1.plot(sat, perm_water, color = 'k', linestyle = '-', marker = 'o')
p2, = ax1.plot(sat, perm_air, color = 'k', linestyle = '-', marker = 'D')
p3, = ax1.plot(x_values1, [x**(3) for x in x_values1], 'k--')
ax1.plot(x_values1, [(1-x)**(3) for x in x_values1], 'k--')
ax1.set_title("strict permeability")
ax1.set_ylabel('permeability')
ax1.set_ylim([0,1])
ax1.set_xlim([0,1])

lgd1 = ax1.legend([p1, p2, p3],
           ["KrWater", "KrAir (strict)", "a = 3"], loc='center left', bbox_to_anchor=(1, 0.5), prop = fontP)

p4, = ax2.plot(sat, diff_water, color = 'k', linestyle = '-', marker = 'o')
p5, = ax2.plot(sat, diff_air, color = 'k', linestyle = '-', marker = 'D')
p6, = ax2.plot(x_values1, [x**(2) for x in x_values1], 'k--')
ax2.plot(x_values1, [(1-x)**(2) for x in x_values1], 'k--')
ax2.set_title("strict diffusivity")
ax2.set_ylabel("diffusivity")
ax2.set_xlabel("saturation")
ax2.set_ylim([0,1])
ax2.set_xlim([0,1])

lgd2 = ax2.legend([p4, p5, p6],
           ["DrWater", "DrAir, (strict)", "a = 3"], loc='center left', bbox_to_anchor=(1, 0.5), prop = fontP)

fig.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1) 
         
fig.show()


# The above code should generate graphs looking something like this:

# In[19]:

i = Image(url = 'http://i.imgur.com/itbMroc.png')
display(i)
g = Image(url = 'http://i.imgur.com/VR8p16K.png')
display(g)


### Completely recreating Gostick's graphs

# The following code more completley recreates Gostick's data.  It checks many more points during saturation to give more data points on the graph, and calculates effective properties in the x, y, and z directions.  There are 4 subplot graphs.  It takes around 15 minutes to finish running.
# 
# If you are running this in IPython Notebook, you will have to restart the kernel so that python does not complain about objects that have the same name as ones already created.

# In[ ]:

import OpenPNM
import matplotlib.pyplot as plt

Lc = 40.5e-6

#setting up network
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

used_inlets = [inlets[x] for x in range(0, len(inlets), 2)]
    
OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=sgl,loglevel=40)
OP_1.setup(invading_fluid = water, defending_fluid = air, inlets = used_inlets,npts=100)
OP_1.run() 


max_inv_seq = max(OP_1.get_throat_data(prop = 'inv_seq'))
    
for x in range(21):
    OP_1.update(seq = max_inv_seq*(x/20))
    
    print('seq = '+str(round(max_inv_seq*(x/20)))+' Seq out of '+str(round(max_inv_seq))+' total sequences')
    
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
            Stokes_alg_single_phase_air.run(loglevel = 40)
            Stokes_alg_single_phase_water.run(loglevel = 40)
            Fickian_alg_single_phase_air.run(loglevel = 40)
            Fickian_alg_single_phase_water.run(loglevel = 40)
            
            Stokes_alg_multi_phase_air.run(loglevel = 40)
            Stokes_alg_multi_phase_water.run(loglevel = 40)
            Fickian_alg_multi_phase_air.run(loglevel = 40)
            Fickian_alg_multi_phase_water.run(loglevel = 40)
            
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
ax3.set_xlim([0,1])


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
ax4.set_xlim([0,1])


lgd2 = ax4.legend([p11, p12, p13, p14, p15, p16, p17, p18, p19, p20],
           ["DrWater,x", "DrWater,y", "DrWater,z",
           "DrAir,x (loose)","DrAir,y (loose)","DrAir,z (loose)",
           "DrAir,x (strict)","DrAir,y (strict)","DrAir,z (strict)", "a = 3"], loc='center left', bbox_to_anchor=(1, 0.5), prop = fontP)

fig.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1) 
         
fig.show()


# With enough patience, the above code should generate a graph looking something like this:

# In[3]:

i = Image(url = 'http://i.imgur.com/sVQiHkM.png')
display(i)


### Discrepancies with Gostick's simulation

# Several things contribute to slight differences between this simulation and that produced by Gostick et al in their 2007 paper.  These include:
# 
# 1) lack of pore size correlation
# 2) lack of late pore filling
