
# coding: utf-8

## Using OpenPNM to recreate Wu's "Determination of oxygen effective diffusivity in porous gas diffusion layer using a three-dimensional pore network model"

# We are going to use OpenPNM to regenerate data found by Wu et al in the paper "Determination of oxygen effective diffusivity in porous gas diffusion layer using a three-dimensional pore network model".  This will serve as a introduction to OpenPNM, and instill an ease with the structure and use of OpenPNM for future simulations.  

### General Structure

# The general procedure we will follow when using OpenPNM to complete this simulation will be:
# 
# 1) create network
# 2) create a geometry object and add methods to calculate gemoetrical properties
# 3) create fluid objects and add methods to calculate fluidic properties
# 4) create physics objects for each fluid and add methods for calculating physical properties for each fluid
# 5) initiate, setup, and run our invasion percolation algorithm
# 6) initiate, setup, and run our Fickian Diffusion Algorithm

### Getting started

# It is assumed that the reader has already downloaded or cloned the most recent version of OpenPNM available through github at https://github.com/PMEAL/OpenPNM.  To make importing OpenPNM as simple as possible, we can add OpenPNM to our pythonpath.  This can be done using the following line of code typed into the command line (specific for MasOSX, and replace /path/to/OpenPNM with the proper path where OpenPNM is saved):  

# In[ ]:

export PYTHONPATH=$PYTHONPATH:/path/to/OpenPNM


# Now that we have added OpenPNM to the python path, we can import it.  We will also need matplotlib.pyplot if we wish to generate graphs of our data in IPython Notebook, and scipy for a few mathematical functions we will use later.

# In[1]:

cd OpenPNM


# In[2]:

import OpenPNM
import matplotlib.pyplot as plt
import scipy as sp


### Generating the pore network and adding a geometry object

# The first step will be to create a pore network we wish to run simulations on.  To regenerate Wu's data, we need to create a pore network that is nxnx2n, where n can be 8, 10, 12, 14, 16, 18, or 20. Wu also specifies that the lattice parameter of the network should be 25e-6.  OpenPNM makes creating a pore network easy.  First, we select the desired network topology (in this case, cubic), then call the generate() method.  Note that divisions[ax,by,cz] supply the number of pores in each of the x, y, and z directions of the lattice.  

# In[3]:

n = 8
Lc = 25e-6

pn = OpenPNM.Network.Cubic(name = 'Wu')

#code to run if we want to set add_boundaries to be False
pn.generate(divisions = [n,n,2*n], add_boundaries= False, lattice_spacing = [Lc])

#code to run if we want to set add_boundaries to be True
#pn.generate(divisions = [n,n,2*n], add_boundaries= True, lattice_spacing = [Lc])


# OpenPNM makes it very easy to visualize the network we have generated through the "Visualization" methods.  We can create vtk files to be viewed using ParaView (downloadable at http://www.paraview.org/download/.  It is suggested that version 3.98 is downloaded instead of 4.1).  If we were able to visualize our pore network model it would appear like this:

# In[4]:

from IPython.display import display
from IPython.display import Image

i = Image(url = 'http://i.imgur.com/ILg7ZdJ.png')
display(i)


# Next, we must create a geometry object so that each pore and throat can be given properties.  For pores, this means diameter and volume, while for throats this means diameter, length, and volume.  We ensure that the geometry will be set for all locations of the network by sending all pore locations as a parameter.
# 
# Afer setting up the geometry object we are ready to add methods to our geometry object for calculating pore and throat dimensional values.  After we add all the methods, we must remember to "regenerate" so that these values are calculated and stored.  Note that if we don't want the logger to show us its progress, we can add 'loglevel = 30' to the parameters of each function.
# 
# The order in which we add each method is important, as the 'sphere' model for calculating pore volumes requires the diameters to be already known.  Similarly, throat lengths and diameters must be set before the throat volumes can be calculated.

# QUESTION: WILL UNIFORM DISTRIBUTION BE INSIDE OPENPNM OR IS THE WAY I HAVE SET THE THROAT DIAMETERS ACCEPTABLE?

# In[5]:

#code to run if boundaries was set to false
geo = OpenPNM.Geometry.GenericGeometry(name = 'wu_geometry', network = pn)
geo.set_locations(pores = pn.pores('all'), throats = 'all')

#code to run if boundaries was set to True
#pn.generate(divisions = [n,n,2*n], add_boundaries= True, lattice_spacing = [Lc], loglevel = 30)
#geo = OpenPNM.Geometry.GenericGeometry(name = 'wu_geometry', network = pn)
#geo.set_locations(pores = pn.pores('internal'), throats = 'all')
#boun = pn.add_geometry(subclass='Boundary',name='boun')
#boun.set_locations(pores=pn.pores('boundary'))

low = .5e-6
high = 9.5e-6

geo.add_method(prop='pore_diameter',model='constant', value = 24e-6)
geo.add_method(prop='pore_volume',model='sphere')
geo.add_method(prop='throat_length',model='straight')
geo.add_method(prop='throat_volume',model='cylinder')

#setting throat diameters to be a uniform distribution
radii = low + sp.rand(pn.num_throats())*(high - low)
pn.set_data(prop = 'diameter', throats = pn.throats(), data = radii*2)

pn.regenerate_geometries()


# Now we can use methods in Tools to return information about our pores and throats.  Note that the printed throat diameters are between 1e-6 and 1.9e-5 (twice our chosen minimum and maximum radii), so we know that the uniform_distribution method is working correctly.  We used the 'straight' model for calculating throat lengths, which simply calculated this length based on the pore diameter and lattice parameter.  Because the lattice parameter is 25e-6 and our pore diameters are 24e-6, our throat lengths should all be 1e-6 which we can check by printing these values.

# In[6]:

throat_diameters = pn.get_throat_data(prop = 'diameter') #if you do not specificy locations, all locations are returned
throat_lengths = pn.get_throat_data(prop = 'length')
print(throat_diameters)
print(throat_lengths)


### Adding fluid objects and methods

# Next, we have to set up our fluids.  We could set up air and water as generic fluids and add methods to each, or we could use the air and water fluids that are already exisiting.  We will use the already existing fluids to make our lives easier.  Again, we will use the regenerate method to make sure that the values for the fluids are calculated and set using the methods chosen (which in this case are preset).  

# In[7]:

air = OpenPNM.Fluids.Air(network = pn, name = 'air')
water = OpenPNM.Fluids.Water(network = pn, name = 'water')
pn.regenerate_fluids()


### Adding Physics objects and methods

# The next step will be to set up physics objects and add the proper methods.  However, Wu does not use the simple bulk_diffusion model we already have to calculate diffusive conductance.  The method we currently have calculates the diffusive conductance accross a conduit (half of one pore, the connecting throat, and half of the next pore) instead of just accross a throat.  We are assuming that Wu calculates the diffusive conductance simply accross a throat.    
# 
# Before we add methods to our physics objects, we should write a method for calculating diffusive conductance that follows Wu's model.  This will not be encorporated into OpenPNM, but can still be used to calculate our diffusive_conductance.  We will call this bulk_diffusion_wu, and it will appear as follows:

# In[8]:

def bulk_diffusion_wu(physics,
                      network,
                      fluid,
                      geometry,
                      propname,
                      diffusivity = 'diffusivity',
                      molar_density = 'molar_density',
                      throat_diameter = 'diameter',
                      throat_length = 'length',
                      pore_diameter = 'diameter',
                      **params):
    r"""
        Calculate the diffusive conductance of throats in network (instead of a
        conduit) based on the areas
        
        Parameters
        ----------
        network : OpenPNM Network Object
        
        fluid : OpenPNM Fluid Object
        The fluid of interest
        
        Notes
        -----
        This function requires that all the necessary fluid properties already be
        calculated.
        
        """
    #ct = fluid.get_data(prop='molar_density',throats='all',mode='interpolate')
    #Interpolate pore values to throats
    DABt = fluid.get_data(prop='diffusivity',throats='all',mode='interpolate')
    #Find g for full throat
    tdia = network.get_throat_data(prop=throat_diameter)
    tlen = network.get_throat_data(prop=throat_length)
    gt = (sp.pi*DABt*tdia**2)/(tlen*4)
    g = gt[geometry.throats()]
    fluid.set_data(prop=propname,throats=geometry.throats(),data=g)


# Now we can set up our physics objects and add the methods with the desired models.  Then we need to "regenerate" to make sure these values are calculated.  Air's diffusive conductance is set a little differently, because we have written our own method that does not exist in OpenPNM.  Note that calling regenerate_physics() will not re-calculate air diffusive conductance.

# In[9]:

phys_water = OpenPNM.Physics.GenericPhysics(network=pn,fluid=water, geometry = geo, name='standard_water_physics')
phys_air = OpenPNM.Physics.GenericPhysics(network=pn,fluid=air, geometry = geo, name='standard_air_physics')

phys_water.add_method(prop='capillary_pressure', model='washburn') #accounts for cylindrical throats
phys_water.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
phys_water.add_method(prop='diffusive_conductance', model='bulk_diffusion', shape = 'circular')
phys_air.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')

bulk_diffusion_wu(physics = phys_air, network = pn, fluid = air, geometry = geo, propname = 'diffusive_conductance')
pn.regenerate_physics()


### Running the invasion percolation algorithm

# Now we are ready to run invasion percolation.  We need to set up our parameters, namely inlets, outlets, and end_condition.  The end_condition parameter, 'breakthrough', specifies that invasion percolation will stop as soon as 1 of the outlet pores is filled.  If we specified 'total' instead, invasion percolation would proceed until all outlet pores are full.
# 
# We will save a vtk file with the saved data so we can view what we've done.  To view this, you can open the file inside ParaView.  

# In[10]:

inlets = pn.get_pore_indices(labels = ['bottom']) #can put in brackets so the whole bottom of the lattice is considered 1 inlet
outlets = pn.get_pore_indices(labels = ['top'])

IP_1 = OpenPNM.Algorithms.InvasionPercolation(network = pn, name = 'OP_1',loglevel=30)
IP_1.setup(invading_fluid = water, defending_fluid = air, inlets = inlets, outlets = outlets, end_condition = 'total')
IP_1.run()
IP_1.update()
#should be uncommented if want to run fickian diffusion on an air-filled lattice
#IP_1.update(IPseq = 0)  

vis = OpenPNM.Visualization.VTK()
vis.write(filename = 'test.vtp', network=pn,fluids=[air,water])


# If all has gone well, we should be able to watch the invasion percolation take place by opening our vtk file in ParaView.  After adding a threshold, we can watch an animation of the invasion that will appear as follows:

# In[11]:

from IPython.display import YouTubeVideo
display(YouTubeVideo('0iSuypRaT7A'))


### Conduit Conductance

# Before we run Fickian diffusion, there is one thing we must take care of.  Throat diffusive conductances need to be reset after invasion percolation is run.  This is because a throat should have a lower air diffusive conductance if water is creating a barrier to air movement.  To deal with this, OpenPNM uses a property called "conduit_conductance".  This updates either diffusive conductance or hydraulic conductance after invasion has taken place (in our case, we want to change diffusive conductance).  This property needs to be added and regenerated after invasion.

# In[12]:

phys_air.add_property(prop='multiphase',model='conduit_conductance',
                  conductance = 'diffusive_conductance', prop_name='conduit_diffusive_conductance',mode='loose')
phys_water.add_property(prop='multiphase',model='conduit_conductance',
                  conductance = 'diffusive_conductance', prop_name='conduit_diffusive_conductance',mode='loose')
phys_air.add_property(prop='multiphase',model='conduit_conductance',
                  conductance = 'hydraulic_conductance', prop_name='conduit_hydraulic_conductance',mode='loose')
phys_water.add_property(prop='multiphase',model='conduit_conductance',
                  conductance = 'hydraulic_conductance', prop_name='conduit_hydraulic_conductance',mode='loose')
pn.regenerate_physics()


### Running the Fickian Diffusion algorithm

# Next, we can set up our Fickian Diffusion algorithm.  Before we run it, we also need to set up our top and bottom boundary.  To continue following Wu's paper, we want to make the top boundary a plane 1/4 from the top, and the bottom boundary 1/4 from the bottom.  This way the Fickian Algorithm is only calculated on the center half.
# 
# If we want our boundaries to be used, we must set the boundary conditions for our algorithm.  Finally, we can setup our algorithm and run it.

# In[13]:

Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loglevel = 30, loggername = 'Fickian', name = 'fickian_alg', network = pn)

A = pn._Nx**2

z_dimension = int(pn.domain_size(dimension = 'height')/Lc) #number of pores in the z direction
quarter_layer = z_dimension/4 #estimates which layer marks 1/4 up the lattice
pore_number = int(quarter_layer*A) #gives the first pore in the layer 1/4 up the lattice

bottom_boundary = list(range(pore_number, pore_number + A))
top_boundary = list(range(pn.num_pores() - pore_number, pn.num_pores() - pore_number +A))

pn.set_pore_info(label='bound1', locations = bottom_boundary)
pn.set_pore_info(label='bound2', locations = top_boundary)

Fickian_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=bottom_boundary)
Fickian_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.2, pores=top_boundary)

Fickian_alg.setup(conductance = 'conduit_diffusive_conductance',fluid=air)
Fickian_alg.run()


# To check that our boundaries have been set properly, we can view the network in ParaView.  The images of the boundaries should appear as below:

# In[14]:

#top boundary
i = Image(url = 'http://i.imgur.com/RknaFjl.png')
display(i)
#bottom boundary
g = Image(url = 'http://i.imgur.com/oxF407s.png')
display(g)


### Calculating effective diffusivity

# The next step will be to use the effective_diffusivity() method to calculate the effective diffusivity of the network after invasion is at completion.  By setting clean as False, we are insuring that the boundaries we have set will be used as the boundaries (otherwise boundaries will be reset at the surfaces and calc_eff_diffusivity will return a tensor that includes the diffusivity in the X, Y, and Z directions).
# 
# When we print normal_diffusivity, we should get something close to zero.  This is because we have run total invasion, and every pore is filled.  This makes it extremely difficult for oxygen to move throughout the system.

# In[15]:

effective_diffusivity = Fickian_alg.calc_eff_diffusivity(clean = False) 
bulk_diffusivity = air.get_pore_data(prop = 'diffusivity')
normal_diffusivity = effective_diffusivity/bulk_diffusivity
print(normal_diffusivity)


### Calculating saturation

# In Wu's first graph, he compares saturation and effective_diffusivity/bulk_diffusivity.  The only value that we are missing is saturation, which we can acquire with the following loop.  We only wish to find the saturation in the section that we used to calculate effective diffusion.  Therefore, we only check pores between the last pore in the bottom boundary and the first pore in the top boundary.  Note that the volume of throats is neglected.
# 
# As a sanity check, if we print saturation here we should get 1.  This is because Invasion Percolation was run with an end_condition of 'total'.

# In[16]:

final_pores = water.get_pore_data('occupancy')*1
pore_volumes = pn.get_pore_data(prop = 'volume')

saturation = sum(final_pores*pore_volumes)/sum(pore_volumes)
print(saturation)


### Generating a graph of saturation versus normalized diffusivity

# Wu's first graph shows normal diffusivity as a function of saturation.  In order to generate this graph, it would be useful to use our above code but include a for loop that updates the invasion to different steps along the way before solving for normal diffusivity.  This way, we can save normal diffusivity at different saturations without generating multiple networks.
# 
# The way the for-loop works:
# as invasion percolation continues, each pore that is invaded is labeled with it's invasion sequence.  When the final pore in invaded, this is the final invasion sequence.  By updating the invasion to different points in this sequence, we can view the network at different steps during the invasion.
# 
# by running this for loop, we are saving different saturation values into x_values and the corresponding normal diffusivity into y_values.  We can later generate a graph with this new data.
# 
# Note that this generates data from only one randomly generated network.  If we want to generate enough data to completely mimic Wu et al, we should run this loop multiple times.  Otherwise, only running one randomly generated network is enough to convince ourselves of the trend.

# In[17]:

max_inv_seq = max(IP_1.get_pore_data(prop = 'IP_inv_seq'))
x_values = []
y_values = []

for x in range(50):
    IP_1.update(IPseq = max_inv_seq*(x/50.0))
    
    phys_air.add_property(prop='multiphase',model='conduit_conductance',
                      conductance = 'diffusive_conductance', prop_name='conduit_diffusive_conductance',mode='loose')
    phys_water.add_property(prop='multiphase',model='conduit_conductance',
                      conductance = 'diffusive_conductance', prop_name='conduit_diffusive_conductance',mode='loose')
    phys_air.add_property(prop='multiphase',model='conduit_conductance',
                      conductance = 'hydraulic_conductance', prop_name='conduit_hydraulic_conductance',mode='loose')
    phys_water.add_property(prop='multiphase',model='conduit_conductance',
                      conductance = 'hydraulic_conductance', prop_name='conduit_hydraulic_conductance',mode='loose')
    pn.regenerate_physics()
    
    Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian', name = 'Fickian', network = pn)
    
    #set labels for top boundary
    #set labels for bottom boundary
    A = pn._Nx**2
    
    z_dimension = int(pn.domain_size(dimension = 'height')/Lc) #number of pores in the z direction
    quarter_layer = z_dimension/4 #estimates which layer marks 1/4 up the lattice
    pore_number = int(quarter_layer*A) #gives the first pore in the layer 1/4 up the lattice
    
    bottom_boundary = list(range(pore_number, pore_number + A))
    top_boundary = list(range(pn.num_pores() - pore_number, pn.num_pores() - pore_number +A))
    
    Fickian_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=bottom_boundary)
    Fickian_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.2, pores=top_boundary)
    
    Fickian_alg.setup(conductance = 'conduit_diffusive_conductance',fluid=air)
    Fickian_alg.run()
     
    effective_diffusivity = Fickian_alg.calc_eff_diffusivity(clean = False) 
    bulk_diffusivity = air.get_pore_data(prop = 'diffusivity')
    
    #calculation of saturation
    final_pores = water.get_pore_data('occupancy')*1
    pore_volumes = pn.get_pore_data(prop = 'volume')

    saturation = sum(final_pores*pore_volumes)/sum(pore_volumes)
    
    x_values.append(saturation)
    y_values.append((effective_diffusivity/bulk_diffusivity)[0])


# the data now has x values saved to x_values and y values saved to y_values.  We've already imported matplotlib.pyplot, so only a few lines of code are needed to generate a graph of the data.  

# In[18]:

plt.plot(x_values, y_values, 'ro')
plt.title('normalized diffusivity versus saturation')
plt.xlabel('saturation')
plt.ylabel('normalized diffusivity')
plt.show()


# The graph generated thus far should look something like this:

# In[19]:

i = Image(url = 'http://i.imgur.com/dlSJGYd.png')
display(i)


### Generating a graph of N (network size parameter) versus f(epsilon)

# Wu's second graph shows how N affects f(epsilon).  f(epsilon) is the same as normalized diffusivity when saturation is 0, therefore we want to always have IP updated to IPseq = 0 (before any pores have been invaded). We want our for loop to change n continuously instead.  Unfortunately, we need n to be set before we begin anything else, so this for loop must be much bigger than the one we were using before.
# 
# We do not need to run invasion percolation, because we want our network to be completely filled with air.  However, in order to have Fickian diffusion run as expected, we must set occupancy values.
# 
# Again, running this loop once will only generate one data point per value for N.  This can ensure us of the trend, but to copy Wu's data we should run this loop more than once.  Note that completion may take a few minutes, as networks that are bigger take longer and longer to process.
# 
# We will use y_2_values to store our normal diffusivity values so that we can save our data from the previous script.

# In[22]:

y_2_values = []

n_values = [8, 10, 12, 14, 16, 18, 20]

for x in range(5):
    for n in n_values:

        Lc = 25e-6

        pn = OpenPNM.Network.Cubic(name = 'Wu')

        #code to run if we want to set add_boundaries to be False
        pn.generate(divisions = [n,n,2*n], add_boundaries= False, lattice_spacing = [Lc])

        geo = OpenPNM.Geometry.GenericGeometry(name = 'wu_geometry', network = pn)
        geo.set_locations(pores = pn.pores('all'), throats = 'all')

        low = .5e-6
        high = 9.5e-6

        geo.add_method(prop='pore_diameter',model='constant', value = 24e-6)
        geo.add_method(prop='pore_volume',model='sphere')

        #setting throat diameters to be a uniform distribution
        radii = low + sp.rand(pn.num_throats())*(high - low)
        pn.set_data(prop = 'diameter', throats = pn.throats(), data = radii*2)

        geo.add_method(prop='throat_length',model='straight')
        geo.add_method(prop='throat_volume',model='cylinder')

        pn.regenerate_geometries()

        #fluids
        air = OpenPNM.Fluids.Air(network = pn, name = 'air')
        water = OpenPNM.Fluids.Water(network = pn, name = 'water')
        pn.regenerate_fluids()

        #physics objects
        phys_water = OpenPNM.Physics.GenericPhysics(network=pn,fluid=water, geometry = geo, name='standard_water_physics')
        phys_air = OpenPNM.Physics.GenericPhysics(network=pn,fluid=air, geometry = geo, name='standard_air_physics')

        phys_water.add_method(prop='capillary_pressure', model='washburn') #accounts for cylindrical throats
        phys_water.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
        phys_water.add_method(prop='diffusive_conductance', model='bulk_diffusion', shape = 'circular')
        phys_air.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')

        bulk_diffusion_wu(physics = phys_air, network = pn, fluid = air, geometry = geo, propname = 'diffusive_conductance')
        pn.regenerate_physics()

        #Invasion percolation
        inlets = pn.get_pore_indices(labels = ['bottom']) #can put in brackets so the whole bottom of the lattice is considered 1 inlet
        outlets = pn.get_pore_indices(labels = ['top'])

        air.set_data(pores = pn.pores(), prop = 'occupancy', data = 1)
        air.set_data(throats = pn.throats(), prop = 'occupancy', data = 1)

        water.set_data(pores = pn.pores(), prop = 'occupancy', data = 0)
        water.set_data(throats = pn.throats(), prop = 'occupancy', data = 0)

        Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian', name = 'Fickian', network = pn)

        #set labels for top boundary
        #set labels for bottom boundary
        A = pn._Nx**2

        z_dimension = int(pn.domain_size(dimension = 'height')/Lc) #number of pores in the z direction
        quarter_layer = z_dimension/4 #estimates which layer marks 1/4 up the lattice
        pore_number = int(quarter_layer*A) #gives the first pore in the layer 1/4 up the lattice

        bottom_boundary = list(range(pore_number, pore_number + A))
        top_boundary = list(range(pn.num_pores() - pore_number, pn.num_pores() - pore_number +A))

        Fickian_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=bottom_boundary)
        Fickian_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.2, pores=top_boundary)

        Fickian_alg.setup(conductance = 'diffusive_conductance',fluid=air)
        Fickian_alg.run()

        effective_diffusivity = Fickian_alg.calc_eff_diffusivity(clean = False) 
        bulk_diffusivity = air.get_pore_data(prop = 'diffusivity')

        y_2_values.append((effective_diffusivity/bulk_diffusivity)[0])


# Now we plot n_values versus y_values, so that we get a graph of N versus f(epsilon).  we can use plt.axis(x_min, x_max, y_min, y_max) to control the axis on the graph, to make it more clear that f(epsilon) should not change with N.  The graph should also make clear that the standard deviation of f(epsilon) decreases as N increases (if you generate more than one data point per n value).  This makes sense, as increasing N also increases the number of throats.
# 
# The plot axis has been adjusted to make it more obvious that F(epsilon) barely changes with N.  

# In[27]:

plt.plot(n_values + n_values + n_values + n_values + n_values, y_2_values, 'ro')
plt.title('F(epsilon) versus N')
plt.xlabel('N')
plt.ylabel('F(epsilon)')
plt.axis(xmin = 6,xmax = 22,ymin= 0,ymax = .2)
plt.show()


# This graph should match the one below: 

# In[29]:

i = Image(url = 'http://i.imgur.com/ktQTci6.png')
display(i)


### Generating a graph of saturation versus g(s)

# Wu's third graph plots saturation versus g(s).  g(s)f(epsilon) = normalized_diffusivity, so g(s) = normalized_diffusivity/f(epsilon).  Because we are not varying the method of generating our network, f(epsilon) will be constant.  Luckily, our second graph calculates this value many times for us.  We should use the average of this value in our calculation of g(s).  Lastly, we can graph g(s) using our x_values from our first graph, and our g(s) values calculated from our y_values and average f(epsilon).    

# In[31]:

#find average value for f(epsilon)
average_f = sum(y_2_values)/len(y_2_values)

#prints graph for g(s) 
g_values = y_values/average_f
 
plt.plot(x_values, g_values, 'ro')
plt.title('g(s) versus saturation')
plt.xlabel('saturation')
plt.ylabel('g(s)')
plt.show()


# The graph generated should mimic the one below.  Note that if the data go above 1 on the g(s) axis this is because we sped up this process by generating less data than Wu.

# In[32]:

i = Image(url = 'http://i.imgur.com/4m2Olmg.png')
display(i)


### Generating a graph of saturation versus alpha

# Finally, we want to print a graph of saturation versus alpha.  This gives us an idea of the equation for g(s) (g(s) = (1-s)^alpha).  The following prints an alpha value for every data point we have gathered.  Wu et al find alpha by using a best fit curve, but this simplification is enough to show that we have alpha values close to those Wu generated.  

# In[33]:

alpha = np.log(1-g_values)*(np.log(x_values))**(-1)
plt.plot(x_values, alpha, 'ro')
plt.title('saturation versus alpha')
plt.xlabel('saturation')
plt.ylabel('alpha')
plt.show()


# The graph generated may look something like this:

# In[34]:

i = Image(url = 'http://i.imgur.com/Sxzvguz.png')
display(i)


### Discrepencies between our data and Wu's data

# Our values aren't quite the same as those generated by Wu et al.  there are several reasons for this-
# 
# 1) the effective diffusivity calculations done by OpenPNM assume a binary system.  Our code finds diffusion of oxygen through nitrogen, while Wu et al find diffusion for oxygen alone.  Being able to completely copy Wu's equations would require much editing to OpenPNM, and it has been decided that maintaining only a binary system way of calculating effective diffusivity is preferable.
# 
# 2) There are some assumptions that had to be made that we cannot be sure matched perfectly with Wu's assumptions.  One example of this is temperature.  We have assumed the oxygen Wu was running through the network was at 273 degrees Kelvin, but we cannot be certain that this is the value he was using.
# 
# 3) There are some slight calculation differences Wu uses because his throats are circular instead of square in cross-section and his pores are spherical instead of cubic.  This mostly effects fluid dynamics calculations including diffusive conductance, hydraulic conductance, and capillary pressure.
