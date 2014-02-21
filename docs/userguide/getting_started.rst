###############################################################################
Getting Started
###############################################################################

===============================================================================
Example: Working with Scripts
===============================================================================
The first thing you must do is import the OpenPNM code so you have access to the functions and methods, so in a blank .py file start by adding the following line at the top:

.. code-block:: python

   import OpenPNM

Next, it's time to generate a network.  This is accomplished by choosing the desired network topology (e.g. cubic), then calling it's `generate` method with the desired parameters:

.. code-block:: python

	pn = OpenPNM.Network.Cubic(name='cubic_1').generate(divisions=[35,35,35],lattice_spacing=[0.0001])

This generates a topological network called `pn` which contains pores at the correct spatial positions and connections between the pores according the desired topology.  

The network does not contain any information about pore and throat sizes at this point.  The next step is to create a geometry object to calculate the desired geometrical properties.  

.. code-block:: python

	geom = OpenPNM.Geometry.GenericGeometry(network=pn,name='stick_and_ball') #instantiate geometry object
	
This freshly instantiated object contains no methods for actual geometry calculations as yet.  A fully functional object is built by adding the desired methods.  For example, the most basic type of geometry is the so-called 'stick and ball' model, where pores are treated as spheres and throats as cylinders.  Furthermore, it is common to assign pore sizes without regard for spatial correlation, but then to assign throat sizes based on the size of the pores it connects.  This is accomplished by choosing the desired models for each property, then adding them to the geometry object.  

.. code-block:: python

	geom.add_method(prop='pore_seed',model='random') #begin adding the desired methods to 'geom'
	geom.add_method(prop='throat_seed',model='neighbor_min')
	geom.add_method(prop='pore_diameter',model='sphere',name='weibull_min',shape=2.5,loc='6e-6',scale=2e-5)
	geom.add_method(prop='throat_diameter',model='cylinder',name='weibull_min',shape=2.5,loc='6e-6',scale=2e-5)
	geom.add_method(prop='pore_volume',model='sphere')
	geom.add_method(prop='throat_volume',model='cylinder')
	geom.add_method(prop='throat_length',model='straight')

OpenPNM ships with many pre-written models available for each property, but adding custom models and even custom properties is designed to be easy.  

At this point the model is now topologically and geometrically complete.  It has pore coordinate, pore and throat sizes and so on.  In order to perform any simulations, however, it is necessary to build fluid objects.  This is done using the same composition technique used to build the geometry.  Fluid objects are instantiated and attached to the network as follows:

.. code-block:: python

	air = OpenPNM.Fluids.GenericFluid(network=pn,name='air')
	water = OpenPNM.Fluids.GenericFluid(network=pn,name='water')
	
Now it is necessary to fill out these two objects with the desired property calculation methods.  For instance, these fluids have a very different viscosity and these must be calculated differently.  

.. code-block:: python

	air.add_method(prop='diffusivity',model='Fuller',MA=0.03199,MB=0.0291,vA=16.3,vB=19.7)
	air.add_method(prop='viscosity',model='Reynolds',uo=0.001,b=0.1)
	air.add_method(prop='molar_density',model='ideal_gas',R=8.314)

	water.add_method(prop='diffusivity',model='constant',value=1e-12)
	water.add_method(prop='viscosity',model='constant',value=0.001)
	water.add_method(prop='molar_density',model='constant',value=44445)
	water.add_method(prop='surface_tension',model='constant',value=0.072)
	water.add_method(prop='contact_angle',model='constant',value=110)
	
The above lines retrieve the requested property estimation method and assign that method to the corresponding property of the fluids.  Any additional arguments or parameter received by these methods are embedded in the function as well.  To determine the surface tension of water now only requires writing `water.surface_tension()`.  Because the model chosen for surface tension was `constant` this method will always return 0.072.  Some of the other models, such as the fuller model of diffusivity, return temperature and pressure dependent values.  If called at this point, `air.diffusivity` will fail because the critical values for the fluid have not been assigned.  This is accomplished using the `set_pore_data` method that is available to the fluid:

.. code-block:: python	

	air.set_pore_data(prop='Pc',data=132.65)
	air.set_pore_data(prop='Tc',data=3.771e6)
	air.set_pore_data(prop='MW',data=0.0291)
	water.set_pore_data(prop='Pc',data=132.65)
	water.set_pore_data(prop='Tc',data=3.771e6)
	water.set_pore_data(prop='MW',data=0.0291)

The above lines add the named properties to the fluid.  Other methods that require such information will now find it when they look for it.  

We are still not ready to perform any experiments, despite the fact that fluids are defined fully built up.  The last step is to define the desired pore scale physics, which defines how the fluid and solid objects interact.  The main example of this is the Washburn equation which predicts the pressure required to push a non-wetting fluid through a capillary of known size.  Defining these pore scale physics models is what differentiates pore network model from each other.  OpenPNM attempts to permit a high degree of extensibility in this area using the same object construction approach used for geometry and fluid above.  Because the physic object defines the interaction of a fluid with the geometry, it is necessary to build one physics object for each fluid:

.. code-block:: python

	phys_water = OpenPNM.Physics.GenericPhysics(network=pn,fluid=water,name='standard_water_physics')
	phys_air = OpenPNM.Physics.GenericPhysics(network=pn,fluid=air,name='standard_air_physics')

As with fluids and geometry objects, the next step is to build-up the bare objects with the desired methods:

.. code-block:: python

	phys_water.add_method(prop='capillary_pressure',model='purcell',r_torioid='1.e-5')
	phys_water.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
	phys_water.add_method(prop='diffusive_conductance',model='bulk_diffusion')

	phys_air.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
	phys_air.add_method(prop='diffusive_conductance',model='bulk_diffusion')

At this point, the system is fully defined and ready for action.  A typical algorithm used in pore network modeling is to use ordinary percolation to simulate drainage of wetting phase by invasion of a nonwetting phase.  An algorithm object must be created as follows:

.. code-block:: python

	OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=pn,name='OP_1')

To perform simulations using this algorithm simply call the `run` command with the desired parameters:

.. code-block:: python
	
	injection_sites = pn.get_pore_indices(subdomain='bottom')
	OP_1.run(invading_fluid='water',defending_fluid='air',inlets=injection_sites,npts=20)
	
The first line in the above block finds all the pores in the network that are labeled 'bottom'.  This labeling step was applied during the network construction.  The list of pores which are to be considered as fluid inlets along with which fluids are the invader and defender are set to the `run` method and the algorithm proceeds.  Upon completion one can view resultant capillary pressure curving using `OP_1.plot_drainage_curve`.

The results of this (and all simulations) are stored locally on the algorithm object.  If these results are desired for use by the rest of the simulation, for subsequent simulations, then they must be explicitly sent out.  This is to prevent unintentional overwriting of results by subsequent algorithms.  This is done using:

.. code-block:: python
	
	OP_a.set_results(Pc=3000)

The above command outputs data called 'occupancy' to the invading fluid object. This data describes which pores and throats are filled by invading and defending fluid at an applied capillary pressure of 3000.  This information can be used by subsequent algorithms.  For instance it is often of interest to determine the gas phase diffusivity through a partially water filled network.  The Fickian diffusion algorithm then would then use this information and set gas diffusion through water filled pores to zero and a relative effective diffusivity value could be found.  

There are many features, details and nuances of this package that have been glossed over in this guide.  The complete documentation describes the OpenPNM framework in detail.  Happy coding.  




