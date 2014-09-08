.. _gostick:

###############################################################################
Example: Regenerating Data from `J.T. Gostick et al. / JPS 173 (2007) 277–290`_
###############################################################################

.. _J.T. Gostick et al. / JPS 173 (2007) 277–290: http://www.sciencedirect.com/science/article/pii/S0378775307009056

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Getting Started
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

In this tutorial, we will regenerate data from J.T. Gostick's 2007 paper `[1]`_. This will both show that OpenPNM can recreate results accurately, and will also show some more specific uses of OpenPNM. While this paper deals with both SGL and Toray GDLs, we will deal only with SGL.

.. _[1]: http://www.sciencedirect.com/science/article/pii/S0378775307009056

There will be a general layout to complete this simulation: 

1. Set up network 
2. Set up geometry and geometrical methods 
3. constrict throat's by a constriction factor 
4. Set up phases and methods 
5. Set up phase physics and methods 
6. Run invasion percolation 
7. Run Stokes and Fickian algorithms 
8. generate effective permeability and effective diffusivity values at different saturations 
9. plot generated data

We first import the OpenPNM code.

.. code-block:: python
    
    import OpenPNM
   
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Setting up Network and Geometry
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

To begin our simulation, we must first generate our SGL network and geometry.  This includes:

1. creating a cubic network object and an SGL10 geometry object
2. sending our geometry object our internal pores
3. calculating values for throat and pore properties for both internal and boundary pores
4. accounting for pores and throats that are too big (making maximum pore size the lattice parameter)

.. code-block:: python

    Lc = 40.5e-6 #Lattice constant used in [1] for SGL 10BA
    #set up network "sgl"
    sgl = OpenPNM.Network.Cubic([26, 26, 10], spacing=Lc, name='sgl')
    sgl.add_boundaries()
    
    #set up geometries, "geo" and "boun"
    Ps = sgl.pores('boundary',mode='difference')
    Ts = sgl.find_neighbor_throats(pores=Ps,mode='intersection',flatten=True)
    geo = OpenPNM.Geometry.SGL10(network=sgl,pores=Ps,throats=Ts,name='geo')
    
    Ps = sgl.pores('boundary')
    Ts = sgl.find_neighbor_throats(pores=Ps,mode='not_intersection')
    boun = OpenPNM.Geometry.Boundary(network=sgl,pores=Ps,throats=Ts,name='boun')
	
Before we move on to setting up our fluid and physics objects, we must constrict throats in the z and y direction by a factor (Gostick et al included this tightening of throats in only these two directions to create realistic anisotropy in the model).  For his SGL simulation, Gostick uses a constriction factor of .95.  Finally, because we have changed values for pore and throat diameters (first by accounting for pores and throats that are too big, and the finally constricting throats in the y and z directions), we must recalculate all pore and throat values relying on these diameters.
	
.. code-block:: python

    #constrict throats in the y and z directions
    throats = sgl.throats('geo')
    connected_pores = sgl.find_connected_pores(throats)
    x1 = [sgl['pore.coords'][pair[0]][0] for pair in connected_pores]
    x2 = [sgl['pore.coords'][pair[1]][0] for pair in connected_pores]
    same_x = [x - y == 0 for x, y in zip(x1,x2)]
    factor = [s*.95 + (not s)*1 for s in same_x]
    throat_diameters = sgl['throat.diameter'][throats]*factor
    
    #remove the regeneration ability of the diameter pore and throat properties
    geo.remove_model('pore.diameter')
    geo.remove_model('throat.diameter')
    boun.remove_model('pore.diameter')
    boun.remove_model('throat.diameter')
    
    #reset aspects relying on pore and throat sizes
    geo.regenerate()
    boun.regenerate()

OpenPNM makes it very easy to visualize the network we have generated through the "Visualization" methods.  We can create vtk files to be viewed using ParaView (downloadable at http://www.paraview.org/download/ ). It is suggested that version 3.98 is downloaded instead of 4.1).  If we visualize our pore network model it would appear like this (the pores have been visualized using boxes- darker boxes are larger.  Because the network is so big, visualization of the throats has been left out for clarity):
	
.. code-block:: python
	
    import OpenPNM.Utilities.IO as io
    io.VTK.save(network=sgl)
	
An example is seen here:

.. image:: http://i.imgur.com/fPZ8lZK.png
	
	
+++++++++++++++++++++++++++++++++
Setting up the Phases and Physics
+++++++++++++++++++++++++++++++++

Now we are ready to set up our phases (water and air) and the physics corresponding to each of these phases. OpenPNM has built in air and water phases, so we can use those. However, Gostick specifies using a water pore contact angle of 100, so we will reset this value after regenerating our fluids.

.. code-block:: python

    #set up phases
    air = OpenPNM.Phases.Air(network = sgl, name = 'air')
    water = OpenPNM.Phases.Water(network = sgl, name = 'water')
    
    #reset pore contact angle
    water['pore.contact_angle'] = 100
    #remove the 
    water.remove_model('pore.contact_angle')

We are now ready to establish physical properties for our fluid objects. To do this, we will: 1) create physics objects associated with our fluids (by using BasePhyics we don't have to add methods for calculating each property because they are already included) 2) use our regenerate_physics() method to calculate these properties

.. code-block:: python

    #create physics objects associated with our phases
    Ps = sgl.pores()
    Ts = sgl.throats()
    phys_water = OpenPNM.Physics.Standard(network=sgl,phase=water,pores=Ps,throats=Ts,dynamic_data=True,name='standard_water_physics')
    phys_air = OpenPNM.Physics.Standard(network=sgl,phase=air,pores=Ps,throats=Ts,dynamic_data=True,name='standard_air_physics')
	
+++++++++++++++++++++++++++++++++
Running Ordinary Percolation, Fickian Diffusion, and Stokes Flow
+++++++++++++++++++++++++++++++++

Gostick uses ordinary percolation to spread water through his GDL before calculating relative permeability and relative diffusivity.  This way, a graph showing the relationship between saturation and relative permeability and between saturation and relative diffusivity can be created.  

To run our ordinary percolation, we will:

1. pick inlet and outlet pores
2. create an Ordinary Percolation algorithm object
3. setup our algorithm object
4. run our algorithm object
5. call update() so that occupancy of pores and throats for each fluid will be set

.. code-block:: python

    inlets = sgl.pores('bottom_boundary')
    used_inlets = [inlets[x] for x in range(0, len(inlets), 2)]
    
    #using every other pore in the bottom and boundary as an inlet
    #prevents extremely small diffusivity and permeability values in the z direction
    used_inlets = [inlets[x] for x in range(0, len(inlets), 2)]
    
    OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=sgl)
    OP_1.run(invading_phase = water, defending_phase = air, inlets = used_inlets,npts=100)

This algorithm performed a start to finish simulation, which fully flooded the network. The 'update_results()' command can be used to update the phase occupancy values throughout the network. 

.. code-block:: python

    #Update the simulation until saturation is at 50%
    OP_1.update_results(sat=0.5) 
	
If we watch a video of the ordinary percolation taking place (which we can do inside paraview), our video should look something like this:

`test animation'_
.. _test animation: http://youtu.be/Fy3bUNTMTUU

The next step will be to calculate effective diffusivity and permeability at different saturations.  Note that we want to run Fickian diffusion and Stokes flow algorithms at different points within our ordinary percolation process.  OpenPNM has a very helpful update() method for updating the occupancy of pores to their values during a specified part of percolation.  During percolation, each pore is given a sequence value showing when in time it was invaded.  We can send update() a sequence parameter, determining when during the percolation we want to update our pore occupancy to.  

The rest of our code will exist within a loop updating our network to different stages of percolation, so that we may view our relative diffusivity and permeability at different points of saturation.

Before we add in the loop aspect, we will walk through the code that will be inside the loop.  

First, we will want to add a physics property that recalculates diffusive and hydraulic conductance in each throat based on occupancy after ordinary percolation has been run.

.. code-block:: python

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

We can finally instatiate, setup, and run our algorithm objects for Stokes flow and Fickian diffusion.  We want to set up 8 different algorithm objects.

1. Stokes flow, single phase air
2. Stokes flow, multi phase air 
3. Stokes flow, single phase water
4. Stokes flow, multi phase water
5. Fickian diffusion, single phase air
6. Fickian diffusion, multi phase air 
7. Fickian diffusion, sing phase water
8. Fickian diffusion, multi phase water

Note that we want the algorithms that are single phase (where only the specified fluid exists in the network) to help us make our permeability and diffusivity values relative.  Any algorithm that is single phase will use the hydraulic or diffusive conductances before we recalculated based on occupancy.  This calls for our conductance parameter to be 'hydraulic_conductance' or 'diffusive_conductance' instead of 'conduit_hydraulic_conductance' or 'conduit_diffusive_conductance'.  

The need for all these different algorithms can be made clearer by the equation relating effective permeability to the absolute permeability and relative permeability: 

:math:`K_{eff, p}(s_p) = K*K_{r, p}(s_p)`

+-------------------------+----------------------------------+
| Key                     | Description                      |
+=========================+==================================+
| :math:`K_{eff, p}(s_p)` | effective permeability of phase  |
|                         | p as a function of saturation    |
+-------------------------+----------------------------------+
| :math:`K`               | absolute permeability (or single |
|                         | phase permeability)              |
+-------------------------+----------------------------------+
| :math:`K_{r, p}(s_p)`   | relative permeability of phase p |
|                         | as a function of saturation      |
+-------------------------+----------------------------------+

Therefore, relative permeability can be found by dividing the effective permeability by the absolute permeability.  Thus the need for a single phase algorithm (absolute permeability) for every multi phase algorithm (effective permeability).

The same goes for relative diffusivity, which has an very similar equation that looks like this:

.. math::

    D_{eff, p}(s_p) = D*D_{r, p}(s_p)

where the same logic applies.

.. code-block:: python
    
    #setting up the 8 StokesFlow and FickianDiffusion algorithms
    Stokes_alg_single_phase_air = OpenPNM.Algorithms.StokesFlow(name='Stokes_alg_single_phase_air',network=sgl,phase=air)
    Stokes_alg_single_phase_water = OpenPNM.Algorithms.StokesFlow(name='Stokes_alg_single_phase_water',network=sgl,phase=water)
    
    Fickian_alg_single_phase_air = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_single_phase_air',network=sgl,phase=air)
    Fickian_alg_single_phase_water = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_single_phase_water',network=sgl,phase=water)
    
    Stokes_alg_multi_phase_air = OpenPNM.Algorithms.StokesFlow(name='Stokes_alg_multi_phase_air',network=sgl,phase=air)
    Stokes_alg_multi_phase_water = OpenPNM.Algorithms.StokesFlow(name='Stokes_alg_multi_phase_water',network=sgl,phase=water)
    
    Fickian_alg_multi_phase_air = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_multi_phase_air',network=sgl,phase=air)
    Fickian_alg_multi_phase_water = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_multi_phase_water',network=sgl,phase=water)

The algorithms are now instantiated, but have not been run yet. In order to run, they need boundary conditions.

.. code-block:: python
    
    #setting boundary conditions
    BC1_pores = sgl.pores(labels='bottom_boundary')
    BC2_pores = sgl.pores(labels='top_boundary')
    
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
    
Now the code can be run. However, we need to be sure that the proper conduit conductance is being employed.

.. code-block:: python

    #run algorithms with proper conduit conductance
    Stokes_alg_single_phase_air.run(conductance = 'hydraulic_conductance')
    Stokes_alg_single_phase_water.run(conductance = 'hydraulic_conductance')
    Fickian_alg_single_phase_air.run(conductance = 'diffusive_conductance')
    Fickian_alg_single_phase_water.run(conductance = 'diffusive_conductance')
    
    Stokes_alg_multi_phase_air.run(conductance = 'conduit_hydraulic_conductance')
    Stokes_alg_multi_phase_water.run(conductance = 'conduit_hydraulic_conductance')
    Fickian_alg_multi_phase_air.run(conductance = 'conduit_diffusive_conductance')
    Fickian_alg_multi_phase_water.run(conductance = 'conduit_diffusive_conductance')

With the algorithms run, each algorithm can calulate it's own effective property.

.. code-block:: python

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

Try printing some of these values out to see how they differ. Remember, that we've just both single and multiphase performed transport simulations in this material.

+++++++++++++++++++++++++++++++++
Running in a large loop to generate graphs
+++++++++++++++++++++++++++++++++



The code at the bottom of this page can be run independantly to generate the a Gostick-like pore network, and it automatically generates the following comparison figure.

.. image:: http://i.imgur.com/eWIM6s2.png

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Discrepancies with Gostick's simulation
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Several things contribute to slight differences between this simulation and that produced by Gostick et al in their 2007 paper.  These include:

1. lack of pore size correlation
2. lack of late pore filling

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
References
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

`[1]`_ J. T. Gostick et al, "Pore network modeling of fibrous gas diffusion layers for polymer electrolyte membrane fuel cells" Journal of Power Sources, vol. 173, issue 1, pp. 277-290, Nov. 2007.

.. code-block:: python

    import OpenPNM
    import matplotlib.pyplot as plt
    
    Lc = 40.5e-6
    
    #1 setting up network
    sgl = OpenPNM.Network.Cubic([26, 26, 10], spacing=Lc, name='SGL10BA')
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
    
    OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=sgl,loglevel=30)
    OP_1.run(invading_phase = water, defending_phase = air, inlets = used_inlets,npts=100)
    
    sat = []
    perm_air = {'0': [], '1': [], '2': []}
    diff_air = {'0': [], '1': [], '2': []}
    perm_water = {'0': [], '1': [], '2': []}
    diff_water = {'0': [], '1': [], '2': []}
    
    max_inv_seq = max(OP_1['throat.inv_seq'])
    
    num_seq = 20
    for x in range(num_seq+1):
        OP_1.update_results(sat = x/num_seq)
    
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
            Stokes_alg_single_phase_air = OpenPNM.Algorithms.StokesFlow(name='Stokes_alg_single_phase_air',network=sgl,phase=air)
            Stokes_alg_single_phase_water = OpenPNM.Algorithms.StokesFlow(name='Stokes_alg_single_phase_water',network=sgl,phase=water)
            
            Fickian_alg_single_phase_air = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_single_phase_air',network=sgl,phase=air)
            Fickian_alg_single_phase_water = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_single_phase_water',network=sgl,phase=water)
            
            Stokes_alg_multi_phase_air = OpenPNM.Algorithms.StokesFlow(name='Stokes_alg_multi_phase_air',network=sgl,phase=air)
            Stokes_alg_multi_phase_water = OpenPNM.Algorithms.StokesFlow(name='Stokes_alg_multi_phase_water',network=sgl,phase=water)
            
            Fickian_alg_multi_phase_air = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_multi_phase_air',network=sgl,phase=air)
            Fickian_alg_multi_phase_water = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_multi_phase_water',network=sgl,phase=water)
            
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

