.. _gostick:

###############################################################################
Tutorial: Regenerating Data from `J.T. Gostick et al. / JPS 173 (2007) 277–290`_
###############################################################################

.. _J.T. Gostick et al. / JPS 173 (2007) 277–290: http://www.sciencedirect.com/science/article/pii/S0378775307009056

===============================================================================
Getting Started
===============================================================================

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

We first import the OpenPNM code and matplotlib.pyplot so that we can plot our graphs at the end.

.. code-block:: python
    
	import OpenPNM
	import matplotlib.pyplot as plt
   
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

	throats = sgl.throats('geo')
	connected_pores = sgl.find_connected_pores(throats)
	x1 = [sgl['pore.coords'][pair[0]][0] for pair in connected_pores]
	x2 = [sgl['pore.coords'][pair[1]][0] for pair in connected_pores]
	same_x = [x - y == 0 for x, y in zip(x1,x2)]
	factor = [s*.95 + (not s)*1 for s in same_x]
	throat_diameters = sgl['throat.diameter'][throats]*factor
	#remove the regeneration ability of the diameter pore and throat properties
	geo.remove_model(models=['pore.diameter','throat.diameter'])
	boun.remove_model(models=['pore.diameter','throat.diameter'])
	#reset aspects relying on pore and throat sizes
	geo.regenerate()
	boun.regenerate()

OpenPNM makes it very easy to visualize the network we have generated through the "Visualization" methods.  We can create vtk files to be viewed using ParaView (downloadable at http://www.paraview.org/download/ ). It is suggested that version 3.98 is downloaded instead of 4.1).  If we visualize our pore network model it would appear like this (the pores have been visualized using boxes- darker boxes are larger.  Because the network is so big, visualization of the throats has been left out for clarity):
	
.. code-block:: python
	
	import OpenPNM.Utilities.IO as io
	io.VTK.save(network=pn,phases=[air,water])
	
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

	OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=sgl,loglevel=30)
	OP_1.run(invading_phase = water, defending_phase = air, inlets = used_inlets,npts=100)
	OP_1.update_results()
	
If we watch a video of the ordinary percolation taking place (which we can do inside paraview), our video should look something like this:

.. youtube:: https://www.youtube.com/watch?feature=player_embedded&v=Fy3bUNTMTUU