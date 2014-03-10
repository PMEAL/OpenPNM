###############################################################################
Getting Started
###############################################################################

===============================================================================
Example: Building a Cubic Network
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
	pn.find_neighbor_pores(pnums=[1])  # [0,2,11,101]
	pn.get_pore_lables(pnums=[1])  # ['all','bottom','left']
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
	geom.add_method(prop='pore_diameter',model='sphere',name='weibull_min',shape=2.5,loc='6e-6',scale=2e-5)
	geom.add_method(prop='throat_diameter',model='cylinder',name='weibull_min',shape=2.5,loc='6e-6',scale=2e-5)
	geom.add_method(prop='pore_volume',model='sphere')
	geom.add_method(prop='throat_volume',model='cylinder')
	geom.add_method(prop='throat_length',model='straight')
	
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
Add Desired Methods to Fluids
-------------------------------------------------------------------------------
	
As with fluids and geometry objects, the next step is to build-up the bare objects with the desired methods:

.. code-block:: python

	phys_water.add_method(prop='capillary_pressure',model='purcell',r_torioid='1.e-5')
	phys_water.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
	phys_water.add_method(prop='diffusive_conductance',model='bulk_diffusion')
	phys_air.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
	phys_air.add_method(prop='diffusive_conductance',model='bulk_diffusion')
	
The final step is to `regenerate()` the object so that the data is actually calculated.  
	
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Run a Drainage Simulation
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

At this point, the system is fully defined and ready to perform some simulations.  A typical algorithm used in pore network modeling is to use ordinary percolation to simulate drainage of wetting phase by invasion of a nonwetting phase.  An Algorithm object is be created as follows:

.. code-block:: python

	OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=pn,name='OP_1')

To perform simulations using this algorithm simply call the `run()` command with the desired parameters:

.. code-block:: python
	
	injection_sites = pn.get_pore_indices(labels='bottom')
	OP_1.run(invading_fluid='water',defending_fluid='air',inlets=injection_sites,npts=20)
	
The first line in the above block finds all the pores in the network that are labeled 'bottom'.  This labeling step was applied during the network construction.  The list of pores which are to be considered as fluid inlets along with which fluids are the invader and defender are set to the `run()` method and the algorithm proceeds.  Upon completion one can view resultant capillary pressure curving using `OP_1.plot_drainage_curve()`.

-------------------------------------------------------------------------------
Sharing Algorithm Results Throughout the Simulation
-------------------------------------------------------------------------------

The results of the above simulation (and all simulations) are stored locally on the algorithm object.  If these results are to be used in other parts of the simulations, then they must be explicitly sent 'out'.  Keeping the results *silo-ed* in this way prevents unintentional overwriting of results by subsequent algorithms.  This allows for multiple simulations of the same type to be run with different conditions and such.  Sending the results of any simulation 'out' is done by with the `update()` command.  Each algorithm :

.. code-block:: python
	
	OP_1.update(Pc=5000)

The above command outputs data called 'occupancy' to the invading fluid object. This data describes which pores and throats are filled by invading and defending fluid at an applied capillary pressure of 5000.  This information can be used by subsequent algorithms.  For instance it is often of interest to determine the gas phase diffusivity through a partially water filled network.  The Fickian diffusion algorithm then would use this information and set gas diffusion through water filled pores to zero and a relative effective diffusivity value could be found.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Run a Diffusion Simulation a Partially Saturated Network
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Calculating the gas phase diffusivity through a water invading porous medium is one of the main applications of pore networks.  The Fickian diffusion algorithm supplied with OpenPNM is setup and called in much the same way as the ordinary percolation algorithm described above.

-------------------------------------------------------------------------------
Prepare the Algorithm
-------------------------------------------------------------------------------

Firstly, an Algorithm object must be instantiated:

.. code-block:: python

	Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(network=pn,name='Fickian_alg')

Each algorithm performs drastically different functions and calculations so each naturally expect quite different arguments.  The Fickian algorithm needs to know what boundary conditions are prevailing.  These can include Dirchelet, various types of Neumann, reaction rates, and so on.  The lines below outline how to setup Dirchelet conditions on two opposing faces.  Note that the process involves first finding the pore finding the indices of pores laying on the 'top' or 'bottom' face of the domain, then applying a the 'Dirichlet' label, and finally applying the boundary value to those locations. 

.. code-block:: python

	BC1 = pn.get_pore_indices(labels=['top'],mode='intersection')
	Fickian_alg.set_pore_info(label='Dirichlet', locations=BC1)
	Fickian_alg.set_pore_data(prop='BCval', data=0.6, locations=BC1)
	BC2 = pn.get_pore_indices(labels=['bottom'],mode='intersection')
	Fickian_alg.set_pore_info(label='Dirichlet', locations=BC2)
	Fickian_alg.set_pore_data(prop='BCval', data=0.2, locations=BC2)
	
Note that this simulation will run on a Network that has been invaded upto 5000 Pa with water due to the OP_1.update(Pc=5000) command used above.  It is a simple matter to change the network saturation be calling this command with a different applied pressure.  

There are many features, details and nuances of this package that have been glossed over in this quickstart guide.  The complete documentation describes the OpenPNM framework in detail.  Happy coding.  




