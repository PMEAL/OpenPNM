.. _tutorial:

###############################################################################
Getting Started
###############################################################################

===============================================================================
Building a Cubic Network
===============================================================================

The first thing you must do is import the OpenPNM code so you have access to the functions and methods, so in a blank *.py* file or at the python command line, start by entering the following line:

.. code-block:: python

   import OpenPNM
   
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Initialize the Network Topology
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Next, it's time to generate a Network.  This is accomplished by choosing the desired network topology (e.g. cubic), then calling it's `generate()` method with the desired parameters:

.. code-block:: python

	pn = OpenPNM.Network.Cubic(name='cubic_1').generate(divisions=[10,10,10],lattice_spacing=[0.0001],add_boundaries=False)

This generates a topological network called *pn* which contains pores at the correct spatial positions and connections between the pores according the desired topology, but without boundary pores.  The network can be queried for certain topological information such as:

.. code-block:: python

	pn.num_pores()  # 1000
	pn.num_throats()  # 2700
	pn.find_neighbor_pores(pores=[1])  # [0,2,11,101]
	pn.get_pore_lables(pores=[1])  # ['all','bottom','left']
	pn.get_pore_indices(labels=['all','bottom','left'],mode='intersection')  # [0,1,2,3,4,5,6,7,8,9]

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Initialize and Build a Geometry Object
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The network does not contain any information about pore and throat sizes at this point.  The next step is to create a geometry object to calculate the desired geometrical properties.  

.. code-block:: python

	geom = OpenPNM.Geometry.GenericGeometry(network=pn,name='stick_and_ball')  # instantiate geometry object
	
-------------------------------------------------------------------------------
Add Desired Methods to Geometry
-------------------------------------------------------------------------------
	
This freshly instantiated object contains no methods for actual geometry calculations as yet.  A fully functional object is built by adding the desired methods.  For example, the most basic type of geometry is the so-called 'stick and ball' model, where pores are treated as spheres and throats as cylinders.  Furthermore, it is common to assign pore sizes without regard for spatial correlation, but then to assign throat sizes based on the size of the pores it connects.  This is accomplished by choosing the desired models for each property, then adding them to the geometry object.  

.. code-block:: python

	geom.add_method(prop='pore_seed',model='random') #begin adding the desired methods to 'geom'
	geom.add_method(prop='throat_seed',model='neighbor_min')
	geom.add_method(prop='pore_diameter',model='sphere',name='weibull_min',shape=2.5,loc=6e-6,scale=2e-5)
	geom.add_method(prop='throat_diameter',model='cylinder',name='weibull_min',shape=2.5,loc=6e-6,scale=2e-5)
	geom.add_method(prop='pore_volume',model='sphere')
	geom.add_method(prop='throat_length',model='straight')
	geom.add_method(prop='throat_volume',model='cylinder')
	
	
Each of the above commands looks into the submodule associated with the `prop` argument, extracts the method corresponding the `model` argument, assigns the specified parameters, and finally attaches the method to the Geometry object.  

Once the object is constructed it is necessary to invoke or run each of the added methods so that they will actually calculate the pore and throat geometry information.  This is done by running the methods, or using the `regenerate()` method which will automatically call all of the methods added through the `add_method()` command:

.. code-block:: python

	geom.pore_seed()
	geom.throat_seed()
	geom.pore_diameter()
	geom.throat_diameter()
	geom.pore_volume()
	geom.throat_volume()
	geom.throat_length()
	geom.regenerate()  # optionally regenerate all methods at once

OpenPNM ships with many pre-written models available for each property, but adding custom models and even custom properties is designed to be easy.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Create Fluids
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

At this point the model is now topologically and geometrically complete.  It has pore coordinates, pore and throat sizes and so on.  In order to perform any simulations, however, it is necessary to build Fluid objects.  This is done using the same composition technique used to build the Geometry.  Fluid objects are instantiated and attached to the Network as follows:

.. code-block:: python

	air = OpenPNM.Fluids.GenericFluid(network=pn,name='air')
	water = OpenPNM.Fluids.GenericFluid(network=pn,name='water')
	
-------------------------------------------------------------------------------
Add Desired Methods to Fluids
-------------------------------------------------------------------------------
	
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
	
The above lines retrieve the requested property estimation method from the submodule indicated by the `prop` argument, and assign that method to the corresponding property of the fluids.  To determine the surface tension of water now only requires writing `water.surface_tension()`.  Because the model chosen for surface tension was `constant` this method will always return 0.072.  Some of the other models, such as the Fuller model of diffusivity, return temperature and pressure dependent values.  If called at this point, `air.diffusivity` will fail because the critical values for the fluid have not been assigned.  This is accomplished using the `set_pore_data` method that is available to the fluid:

.. code-block:: python	

	air.set_pore_data(prop='Pc',data=132.65)
	air.set_pore_data(prop='Tc',data=3.771e6)
	air.set_pore_data(prop='MW',data=0.0291)
	water.set_pore_data(prop='Pc',data=132.65)
	water.set_pore_data(prop='Tc',data=3.771e6)
	water.set_pore_data(prop='MW',data=0.0291)

The above lines add the named properties to the fluid.  Other methods that require such information will now find it when they look for it.  

Like the Geometry object, it is necessary to actually run each of the added methods for the data to be generated.  This can also be accomplished with the `regenerate()` command.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Create Pore Scale Physics Objects
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

We are still not ready to perform any experiments, despite the fact that fluids are defined fully built up.  The last step is to define the desired pore scale physics, which defines how the fluid and solid objects interact.  A classic example of this is the Washburn equation which predicts the pressure required to push a non-wetting fluid through a capillary of known size.  OpenPNM attempts to permit a high degree of extensibility by using the same object construction approach used for Geometry and Fluid above.  Because the Physics object defines the interaction of a Fluid with the Geometry, it is necessary to build one physics object for each Fluid (and Geometry).  

.. code-block:: python

	phys_water = OpenPNM.Physics.GenericPhysics(network=pn,fluid=water,name='standard_water_physics')
	phys_air = OpenPNM.Physics.GenericPhysics(network=pn,fluid=air,name='standard_air_physics')

-------------------------------------------------------------------------------
Add Desired Methods to Physics Objects
-------------------------------------------------------------------------------
	
As with fluids and geometry objects, the next step is to build-up the bare objects with the desired methods:

.. code-block:: python

	phys_water.add_method(prop='capillary_pressure',model='purcell',r_toroid=1.e-5)
	phys_water.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
	phys_water.add_method(prop='diffusive_conductance',model='bulk_diffusion')
	phys_air.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
	phys_air.add_method(prop='diffusive_conductance',model='bulk_diffusion')
	
The final step is to ``regenerate()`` the object so that the data is actually calculated.  
	
 






