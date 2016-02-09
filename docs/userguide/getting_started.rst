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
	mgr = OpenPNM.Base.Workspace()

The *Workspace* object provides high-level oversight to all the simulations existing in memory at any given time.  Its main purpose is saving, loading and exporting data to files.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Initialize the Network Topology
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Next, it's time to generate a *Network*.  This is accomplished by choosing the desired network topology (e.g. cubic), then calling its respective method in OpenPNM with the desired parameters:

.. code-block:: python

	pn = OpenPNM.Network.Cubic(name='net', shape=[10,10,10])

This generates a topological network called ``pn`` which contains pores at the correct spatial positions and connections between the pores according the desired topology, but without boundary pores.  The network can be queried for certain topological information such as:

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

Networks are the chief objects in OpenPNM.  All subsequent objects are subjects of the Network.  Two Networks may exist in memory at the same time, but they are completely independent of each other.  In fact, because of thier chief status, the *Network* object define a simulation, so the term 'Network' and 'simulation' are sometimes used interchangably. The *Workspace* object mentioned above keeps track of all the different *Networks* in memory.

===============================================================================
Initialize and Build a Geometry Object
===============================================================================

The *Network* does not contain any information about pore and throat sizes at this point.  The next step is to create a *Geometry* object to calculate the desired geometrical properties.

.. code-block:: python

	geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)

This statement contains three arguments.  The ``network`` tells the *Geometry* object which *Network* it is associated with.  ``pores`` and ``throats`` indicate which locations in the *Network* this *Geometry* object will apply to.  It's possible to have multiple *Geometry* Objects each applying to different regions of the domain, to create heterogeneous materials for instance.  *Geometry* objects cannot overlap, however, since it does not make sense for a single pore to have multiple sizes.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Add Desired Properties to Geometry
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This freshly instantiated *Geometry* object contains no geometric properties as yet.  For this example we will make a basic 'stick-and-ball' geometry, which as the name suggests treats the pores as spheres and the throats as cylinders.

-------------------------------------------------------------------------------
Direct Assignment of Static Values
-------------------------------------------------------------------------------

Before diving into the 'pore-scale model' feature, let's assign a random seed to each pore.  These seed values will be subsequently used in the calculation of pore diameters using a statistical distribution.

.. code-block:: python

    geom['pore.seed'] = sp.rand([pn.Np,])

This creates a Numpy array of random numbers that is *Np* long, meaning each pore is assigned a unique random number. This is one of the the most basic ways to assign values to the Geomtry object.  The limitation of this approach is that the values are now completley static and can only be updated by re-assigning random values.  In some cases it may be of interest to have values *regenerate* upon request and this functionality is provided by the *models* feature to be described next.

-------------------------------------------------------------------------------
Assigning Pore Scale Models to Calculate Properties
-------------------------------------------------------------------------------

OpenPNM includes an array of prewritten pore-scale models which are found in the *models* folder under each submodule.  To access these models, the first step is to load the Geometry model library into a convenient namespace.

.. code-block:: python

	import OpenPNM.Geometry.models as gm

The 'behind-the-scenes' behavior that occurs when adding a pore scale model to an object is outlined in it's own  :ref:`documentation page<models>`.  For the purpose of this guide these details will be skipped.  To add a model, you can either use to the ```object.models.add`` or ``object.add_models`` command.  For instance, OpenPNM comes with a model for assigning random values to pores, instead of the direct assignment above:

.. code-block:: python

    geom.add_model(propname='pore.seed', model=gm.pore_misc.random)

The above line generates an *Np* long list of random numbers and insert them into the ``geom['pore.seed']`` exactly we did previously.  The difference is that when we call :code:`geom.regenerate()` the random numbers will be regenerated...and so will all the other values in ``geom`` that are calculated by a pore scale model!  This mechanism enables the changes in one property to cascade to all other relevant properties.

Each pore scale model takes different arguments.  In the code block below, a Weibull distribution is assigned to the pore diameters, which will use the ``'pore.seed'`` values, the throat diameter is taken as the minimum of its two neighbors, and other geoemtric properties are calculated in the expected way.

.. code-block:: python

	geom.add_model(propname='pore.diameter',
                 model=gm.pore_diameter.sphere,
                 psd_name='weibull_min',
                 psd_shape=2.77,
                 psd_loc=6.9e-7,
                 psd_scale=9.8e-6,
                 psd_offset=10e-6)
  geom.add_model(propname='throat.diameter',
                 model=gm.throat_misc.neighbor,
                 pore_prop='pore.diameter',
                 mode='min')
  geom.add_model(propname='pore.volume', model=gm.pore_volume.sphere)
  geom.add_model(propname='pore.area', model=gm.pore_area.spherical)
  geom.add_model(propname='throat.length', model=gm.throat_length.straight)
  geom.add_model(propname='throat.volume', model=gm.throat_volume.cylinder)
  geom.add_model(propname='throat.area', model=gm.throat_area.cylinder)

At this point, ``geom`` has been fully populated with the necessary geometric properties.  You can view these by typing ``print(geom)`` at the command line.

===============================================================================
Create Phases
===============================================================================

The simulation is now topologically and geometrically complete.  It has pore coordinates, pore and throat sizes and so on.  In order to perform any simulations, however, it is necessary to build *Phase* objects that represent the fluids in the simulations.  This is done using the same composition technique used to build the *Geometry*.  Phases objects are instantiated as follows:

.. code-block:: python

	air = OpenPNM.Phases.GenericPhase(network=pn, name='air')
	water = OpenPNM.Phases.GenericPhase(network=pn, name='water')

Again, note ``pn`` is passed as an argument because this *Phase* must know to which *Network* it belongs.  Also, note that ``pores`` and ``throats`` are NOT specified; this is because *Phases* are assumed to exist everywhere in the domain.  For multiphase immiscible flow the presence or absence of a *Phase* in given locations is tracked using a ``'pore.occupancy'`` array.

.. note:: **Naming Objects**

	The above two lines also include a ``name`` argument.  All objects in OpenPNM can be named in this way if desired, however, if no name is given one will be generated.  The point of the name is to allow easy identification of an object at the command line, using the ``name`` attribute (``air.name``).  Objects can be renamed, so if you wish to override a default name simply use ``air.name`` = 'air'.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Add Desired Methods to Phases
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Now it is necessary to fill these two *Phase* objects with the desired thermophysical properties.  For instance, they may have very different viscosity and these must be calculated differently. It is possible to simply 'hard code' static property values, as follows:

.. code-block:: python

	water['pore.temperature'] = 298.0
	water['pore.diffusivity'] = 1e-12
	water['pore.viscosity'] = 0.001
	water['pore.molar_density'] = 44445.0
	water['pore.contact_angle'] = 110.0
	water['pore.surface_tension'] = 0.072

It should be reiterated here that these static property values are not updated when other properties change.  For instance, if the temperature of the simulation is changed to 353 K from 298 K, the viscosity must also change.  Using static values for properties means that viscosity must be recalculated and re-assigned manually.  The 'pore-scale model' approach addresses this.

.. note:: **Scalar to Vector Conversion During Assignment**

	The above block illustrates a feature of OpenPNM that is worth pointing out now.  All pores need to have a diffusivity value associated with them; however, we often want to assign the same value to every pore.  If you assign a scalar value to any property in OpenPNM it will automatically be converted to a vector of the appropriate length (either *Np* or *Nt* long).  This is explained in more detail :ref:`here<inner_workings>`.

To use the available thermophysical property models that are included with OpenPNM, import the *Phase* models library:

.. code-block:: python

	from OpenPNM.Phases import models as fm

The code block below illustrate how to define a *Phase* object to represent Air using 'pore-scale models'. Some of the models require various input parameters.  For instance, consider the Fuller model, which requires the molecular mass and diffusion volume of the species in the mixture.  More importantly, the Fuller model also includes temperature, meaning that if temperature of the phase changes, then the model can be re-run to regenerate the diffusivity at the new temperature.  The Fuller model code assumes that the temperature for the *Phase* can be found in ``'pore.temperature'``.  It's possible to customize these default property names as outlined :ref:`here<customizing>`.

.. code-block:: python

  air.add_model(propname='pore.diffusivity',
                model=fm.diffusivity.fuller,
                MA=0.03199,
                MB=0.0291,
                vA=16.3,
                vB=19.7)
  air.add_model(propname='pore.viscosity',
                model=fm.viscosity.reynolds,
                uo=0.001,
                b=0.1)
  air.add_model(propname='pore.molar_density',
                model=fm.molar_density.ideal_gas,
                R=8.314)

===============================================================================
Create Pore Scale Physics Objects
===============================================================================

We are still not ready to perform any simulations.  The last step is to define the desired pore scale physics, which defines how the phase and geometrical properties interact.  A classic example of this is the Washburn equation which predicts the capillary pressure required to push a non-wetting fluid through a capillary of known size.  Because the *Physics* object defines the interaction of a *Phase* with the *Geometry*, it is necessary to build one *Physics* object for each intersection between *Geometry* and *Phase* objects:

.. code-block:: python

	phys_water = OpenPNM.Physics.GenericPhysics(network=pn,
	                                            phase=water,
                                              geometry=geom)
	phys_air = OpenPNM.Physics.GenericPhysics(network=pn,
	                                          phase=air,
                                            geometry=geom)

*Physics objects do not require the specification of which ``pores`` and ``throats`` they apply.  This assignment is implied by the passing of a ``geometry`` object, which has already been assigned to specific locations.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Add Desired Methods to Physics Objects
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

As with *Phase* and *Geometry* objects, the next steps are to load the model library from the *Physics* submodule and then to build-up the bare objects with the desired models:

.. code-block:: python

	from OpenPNM.Physics import models as pm

	phys_water.add_model(propname='throat.capillary_pressure',
	                     model=pm.capillary_pressure.washburn)
	phys_water.add_model(propname='throat.hydraulic_conductance',
	                     model=pm.hydraulic_conductance.hagen_poiseuille)
	phys_air.add_model(propname='throat.diffusive_conductance',
	                   model=pm.diffusive_conductance.bulk_diffusion)
	phys_air.add_model(propname='throat.hydraulic_conductance',
	                   model=pm.hydraulic_conductance.hagen_poiseuille)


===============================================================================
Run some simulations
===============================================================================

Finally, it is now possible to run some simulations.  The code below estimates the effective diffusivity through the network by applying a concentration gradient across and calculating the flux.  This starts by creating a FickianDiffusion *Algorithm*, which is pre-defined in OpenPNM:

.. code-block:: python

	alg = OpenPNM.Algorithms.FickianDiffusion(network=pn,phase=air)

Next the boundary conditions are applied using the ``set_boundary_conditions`` method.  In this case the boundary conditions are applied to the ``'left'`` and ``'right'`` of the cubic domain.

.. code-block:: python

	# Assign Dirichlet boundary conditions to top and bottom surface pores
	BC1_pores = pn.pores('right')
	alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
	BC2_pores = pn.pores('left')
	alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.4, pores=BC2_pores)

.. note:: **Pore and Throat Labels**

	Note how the ``pores`` method was used to extract pore numbers based on the labels ``'left'`` and ``'right'``.  It's possible to add your own labels to the simulations to allow quick access to special sets of pores.  This is outlined :ref:`here<inner_workings>`.

To actually run the algorithm use the ``run`` method.  This builds the coefficient matrix from the existing values of diffusive conductance, and inverts the matrix to solve for concentration in each pores.

.. code-block:: python

	# Use desired diffusive_conductance in the diffusion calculation (conductance for the dry network or water-filled network)
	alg.run()
	alg.return_results()
	# Calculate the macroscopic effective diffusivity through this Network
	Deff = alg.calc_eff_diffusivity()

===============================================================================
Visualise the Results
===============================================================================
We can now visualise our network and simulation results.  OpenPNM does not support native visualization, so data must be exported to a file for exploration in another program such as any of the several VTK front ends (i.e. Paraview).

.. code-block:: python

	mgr.export(pn)

This creates a *net.vtp* file in the active directory, which can be loaded from ParaView. For a quick tutorial on the use of Paraview with OpenPNM data, see :ref:`Using Paraview<paraview_example>`.

To save an incomplete simulation for later work, the **Workspace** object can be used to save the entire workspace (i.e. all simulations) using ``mgr.save()``, or just the simulation of interest using ``mgr.save_simulation(pn)``.
