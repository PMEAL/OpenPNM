.. _getting_started:

###############################################################################
A Quick Guide to Getting Started with OpenPNM
###############################################################################

===============================================================================
Building a Cubic Network
===============================================================================

The first thing you must do is import the OpenPNM code so you have access to the functions and methods, so in a blank *.py* file or at the python command line, start by entering the following line:

.. code-block:: python

	import OpenPNM
	ctrl = OpenPNM.Base.Controller()

The **Controller** object provides high-level oversight to all the simulations existing in memory at any given time.  Its main purpose is saving, loading and export data to files.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Initialize the Network Topology
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Next, it's time to generate a Network.  This is accomplished by choosing the desired network topology (e.g. cubic), then calling its respective method in OpenPNM with the desired parameters:

.. code-block:: python

	pn = OpenPNM.Network.Cubic(name='net',shape=[10,10,10])

This generates a topological network called *pn* which contains pores at the correct spatial positions and connections between the pores according the desired topology, but without boundary pores.  The network can be queried for certain topological information such as:

.. code-block:: python

	pn.num_pores()  # 1000
	pn.num_throats()  # 2700
	pn.find_neighbor_pores(pores=[1])  # [0,2,11,101]
	pn.labels(pores=[1])  # ['all','bottom','left']
	pn.pores(labels = 'bottom')

This data may also be stored in a variable:

.. code-block:: python

	Ps = pn.pores()
	Ts = pn.throats()

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Initialize and Build a Geometry Object
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The network does not contain any information about pore and throat sizes at this point.  The next step is to create a Geometry object to calculate the desired geometrical properties.

.. code-block:: python

	geom = OpenPNM.Geometry.GenericGeometry(network=pn,pores=Ps,throats=Ts)  # instantiate geometry object

-------------------------------------------------------------------------------
Add Desired Methods to Geometry
-------------------------------------------------------------------------------

This freshly instantiated object contains no methods for actual geometry calculations as yet.  A fully functional object is built by adding the desired methods.  For example, the most basic type of geometry is the so-called 'stick and ball' model, where pores are treated as spheres and throats as cylinders.  Furthermore, it is common to assign pore sizes without regard for spatial correlation, but then to assign throat sizes based on the size of the pores it connects.  This is accomplished by choosing the desired models for each property, then adding them to the geometry object.

The first step is to load the Geometry model library.

.. code-block:: python

	import OpenPNM.Geometry.models as gm

Then, the different geometry models are added one by one to the object geom.

.. code-block:: python

    geom.add_model(propname='pore.seed',model=gm.pore_misc.random) #begin adding the desired methods to 'geom'
    geom.add_model(propname='pore.diameter',
                   model=gm.pore_diameter.sphere,
                   psd_name='weibull_min',
                   psd_shape=2.77,
                   psd_loc=6.9e-7,
                   psd_scale=9.8e-6,
                   psd_offset=10e-6)
    geom.add_model(propname='throat.diameter',model=gm.throat_misc.neighbor,pore_prop='pore.diameter',mode='min')
    geom.add_model(propname='pore.volume',model=gm.pore_volume.sphere)
    geom.add_model(propname='pore.area',model=gm.pore_area.spherical)
    geom.add_model(propname='throat.length',model=gm.throat_length.straight)
    geom.add_model(propname='throat.volume',model=gm.throat_volume.cylinder)
    geom.add_model(propname='throat.area',model=gm.throat_area.cylinder)

Each of the above commands extracts the model, assigns the specified parameters, and attaches the model to the Geometry object.

OpenPNM ships with many pre-written models available for each property, but adding custom models and even custom properties is designed to be easy.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Create Phases
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

At this point the model is now topologically and geometrically complete.  It has pore coordinates, pore and throat sizes and so on.  In order to perform any simulations, however, it is necessary to build Phases objects that e.g. represent the fluids in the simulations.  This is done using the same composition technique used to build the Geometry.  Phases objects are instantiated and attached to the Network as follows:

.. code-block:: python

	air = OpenPNM.Phases.GenericPhase(network=pn,name='air')
	water = OpenPNM.Phases.GenericPhase(network=pn,name='water')

-------------------------------------------------------------------------------
Add Desired Methods to Phases
-------------------------------------------------------------------------------

Now it is necessary to fill out these two objects with the desired property calculation model.  For instance, these phases have a very different viscosity and these must be calculated differently.
As for the geometric object, the phase models need to be load first:

.. code-block:: python

	from OpenPNM.Phases import models as fm

Then, water and air properties are then defined by the code below. Note that some of the models, such as the Fuller model of diffusivity, needs input parameters as molar masses. These inputs are simply state in the add_model method.

.. code-block:: python

    air.add_model(propname='pore.diffusivity',model=fm.diffusivity.fuller,MA=0.03199,MB=0.0291,vA=16.3,vB=19.7)
    air.add_model(propname='pore.viscosity',model=fm.viscosity.reynolds,uo=0.001,b=0.1)
    air.add_model(propname='pore.molar_density',model=fm.molar_density.ideal_gas,R=8.314)
    water['pore.diffusivity'] = 1e-12
    water['pore.viscosity'] = 0.001
    water['pore.molar_density'] = 44445.0
    water['pore.contact_angle'] = 110.0
    water['pore.surface_tension'] = 0.072


The first above lines retrieve the requested property estimation model from the submodule indicated by the `propname` argument, and assign that method to the corresponding property of the phases on each pore location.  The last five lines set a constant value, by placing it directly into a new dictionary entry.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Create Pore Scale Physics Objects
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

We are still not ready to perform any simulations.  The last step is to define the desired pore scale physics, which defines how the phase and geometrical properties interact.  A classic example of this is the Washburn equation which predicts the pressure required to push a non-wetting fluid through a capillary of known size.  Because the Physics object defines the interaction of a Phase with the Geometry, it is necessary to build one physics object for each intersection between Geometry and Phase objects:

.. code-block:: python

	phys_water = OpenPNM.Physics.GenericPhysics(network=pn,phase=water,geometry=geom)
	phys_air = OpenPNM.Physics.GenericPhysics(network=pn,phase=air,geometry=geom)

-------------------------------------------------------------------------------
Add Desired Methods to Physics Objects
-------------------------------------------------------------------------------

As with phases and geometry objects, the next steps are first to load the model library and to build-up the bare objects with the desired models:

.. code-block:: python

	from OpenPNM.Physics import models as pm

	phys_water.add_model(propname='throat.capillary_pressure',
	                     model=pm.capillary_pressure.purcell,
	                     r_toroid=1.e-5)
	phys_water.add_model(propname='throat.hydraulic_conductance',
	                     model=pm.hydraulic_conductance.hagen_poiseuille)
	phys_air.add_model(propname='throat.diffusive_conductance',
	                   model=pm.diffusive_conductance.bulk_diffusion)
	phys_air.add_model(propname='throat.hydraulic_conductance',
	                   model=pm.hydraulic_conductance.hagen_poiseuille)
-------------------------------------------------------------------------------
Run some simulations
-------------------------------------------------------------------------------

.. code-block:: python

	alg = OpenPNM.Algorithms.FickianDiffusion(network=pn,phase=air)
	# Assign Dirichlet boundary conditions to top and bottom surface pores
	BC1_pores = pn.pores('right')
	alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
	BC2_pores = pn.pores('left')
	alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.4, pores=BC2_pores)
	# Use desired diffusive_conductance in the diffusion calculation (conductance for the dry network or water-filled network)
	alg.run()
	alg.return_results()
	# Calculate the macroscopic effective diffusivity through this Network
	Deff = alg.calc_eff_diffusivity()

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Visualise the Results
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
We can now visualise our network and simulation results.  OpenPNM does not (yet) support native visualization, so data must be exported to a file for exploration in another program such as any of the several VTK front ends (I.e. Paraview).

.. code-block:: python

	ctrl.export(pn)

This creates a *net.vtp* file in the active directory, which can be loaded from ParaView. For a quick tutorial on the use of Paraview with OpenPNM data, see :ref:`Using Paraview<paraview_example>`.

To save an incomplete simulation for later work, the **Controller** object can be used to save the entire workspace (i.e. all simulations) using ``ctrl.save()``, or just the simulation of interest using ``ctrl.save_simulation(pn)``.
