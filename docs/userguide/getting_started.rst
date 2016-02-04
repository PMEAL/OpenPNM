.. _getting_started:

###############################################################################
A Quick Guide to Getting Started with OpenPNM
###############################################################################

===============================================================================
Building a Cubic Network
===============================================================================

Start by generating a *Network*.  This is accomplished by choosing the desired network topology (e.g. cubic), then calling its respective method in OpenPNM with the desired parameters:

.. code-block:: python

	pn = OpenPNM.Network.Cubic(shape=[10, 10, 10], spacing=0.0001)

This generates a topological network and stores it in variable ``pn``.  This network contains pores at the correct spatial positions and connections between the pores according the specified topology (but without boundary pores).  The ``shape`` argument specifies the number of pores in the [X, Y, Z] directions of the cube.  Networks in OpenPNM are alway 3D dimensional, meaning that a 2D or 'flat' network is still 1 layer of pores 'thick' so [X, Y, Z] = [20, 10, 1].  The ``spacing`` argument controls the center-to-center distance between pores.  Although OpenPNM does not currently have a dimensional units system, we *strongly* recommend using SI throughout.

The network can be queried for a variety of common topological:

>>> pn.num_pores()  # 1000
>>> pn.num_throats()  # 2700
>>> pn.find_neighbor_pores(pores=[1])  # [0, 2, 11, 101]
>>> pn.labels(pores=[1])  # ['all', 'bottom', 'left']
>>> pn.pores(labels='bottom')

The data returned from these queries may also be stored in a variable for convenience:

>>> Ps = pn.pores()
>>> Ts = pn.throats()

===============================================================================
Initialize and Build a Geometry Object
===============================================================================

The *Network* does not contain any information about pore and throat sizes at this point.  The next step, then, is to create a *Geometry* object to calculate the desired geometrical properties.

>>> geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)

This statement contains three arguments: ``network`` tells the *Geometry* object which *Network* it is associated with.  ``pores`` and ``throats`` indicate which locations in the *Network* where this *Geometry* object will apply.

.. note::

	OpenPNM was designed to allow multiple *Geometry* objects, with each applying to different regions of the *Network*.  This enables modeling of heterogeneous materials with much different geometrical properties in different regions; this is why the ``pores`` and ``throats`` arguments are required.  In this tutorial ``geom`` applies everywhere which is a common scenario.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Add Desired Properties to Geometry
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This freshly instantiated *Geometry* object contains no geometric properties as yet because we chose to use the *GenericGeometry* class.  There are several other 'pre-written' classes available in *Geometry* module, but for this example we will explicity make a 'stick-and-ball' geometry, which as the name suggests treats the pores as spheres and the throats as cylinders.

-------------------------------------------------------------------------------
Direct Assignment of Static Values
-------------------------------------------------------------------------------

Let's start by assiging diameters to each pore from a random distribution, spanning 10 um to 100 um.  The upper limit is arise because the ``spacing`` of the *Network* was set to 100 [um], so pore diameters exceeding 100 um might overlap with neighbors.  The lower limit is merely to avoid vanishingly small pores.

>>> geom['pore.diameter'] = 0.00001 + sp.rand(pn.Np)*0.00099

This creates a ND-array of random numbers (between 0.00001 and 0.0001) that is *Np* long, meaning each pore is assigned a unique random number.

For throat diameter, we want them to always be smaller than the two pores which it connects to maintain physical consistency. This requires explaining how OpenPNM stores network topology.

>>> P12 = pn['throat.conns']  # An Nt x 2 list of pores on the end of each throat
>>> D12 = geom['pore.diameter'][P12]  # An Nt x 2 list of pore diameters
>>> Dt = sp.amin(D12, axis=1)  # An Nt x 1 list of the smaller pore from each pair
>>> geom['throat.diameter'] = Dt

Let's disect the above lines.  Firstly, P12 is a direct copy of the Network's \'throat.conns\' array, which contains the indices of the pore pair connected by each throat.  Next, this *Nt-by-2* array is used to index into the \'pore.diameter'\ array, resulting in another *Nt-by-2* array containing the diameters of the pores connected by each throat.  Finally, the Scipy function ``amin`` is used to find the minimum diameter of each pore pair by specifying the ``axis`` keyword as 1, and the resulting *Nt-by-1* array is assigned to ``geom['throat.diameter']``.

Finally, we must specify the remaining geometrical properties of the pores and throats. Since we're creating a 'stick-and-ball' geometry, the sizes are calculated from the geometrical equations for spheres and cylinders as follows:

>>> Rp = geom['pore.diameter']/2
>>> geom['pore.volume'] = (4/3)*3.14159*(Rp)**3
>>> geom['throat.length'] = ??
>>> Rt = geom['throat.diameter']/2
>>> Lt = geom['throat.length']
>>> geom['throat.volume'] = 3.14159*(R)**2*L

The basic geometrical properties of the network are now defined.

===============================================================================
Create Phases
===============================================================================

The simulation is now topologically and geometrically complete.  It has pore coordinates, pore and throat sizes and so on.  In order to perform any simulations, however, it is necessary to build *Phase* objects that represent the fluids in the simulations.  This is done using the same composition technique used to build the *Geometry*.  Phases objects are instantiated as follows:

>>> air = OpenPNM.Phases.GenericPhase(network=pn, name='air')
>>> water = OpenPNM.Phases.GenericPhase(network=pn, name='water')

Again, note ``pn`` is passed as an argument because this *Phase* must know to which *Network* it belongs.  Also, note that ``pores`` and ``throats`` are NOT specified; this is because *Phases* are assumed to exist everywhere in the domain.  For multiphase immiscible flow the presence or absence of a *Phase* in given locations is tracked using a ``'pore.occupancy'`` array.

.. note:: **Naming Objects**

	The above two lines also include a ``name`` argument.  All objects in OpenPNM can be named in this way if desired, however, if no name is given one will be generated.  The point of the name is to allow easy identification of an object at the command line, using the ``name`` attribute (``air.name``).  Objects can be renamed, so if you wish to override a default name simply use ``air.name`` = 'air'.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Add Desired Methods to Phases
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Now it is necessary to fill these two *Phase* objects with the desired thermophysical properties.  For instance, they may have very different viscosity and these must be calculated differently. It is possible to simply 'hard code' static property values, as follows:

>>> water['pore.temperature'] = 298.0
>>> water['pore.diffusivity'] = 1e-12
>>> water['pore.viscosity'] = 0.001
>>> water['pore.molar_density'] = 44445.0
>>> water['pore.contact_angle'] = 110.0
>>> water['pore.surface_tension'] = 0.072

It should be reiterated here that these static property values are not updated when other properties change.  For instance, if the temperature of the simulation is changed to 353 K from 298 K, the viscosity must also change.  Using static values for properties means that viscosity must be recalculated and re-assigned manually.  The 'pore-scale model' approach addresses this.

.. note:: **Scalar to Vector Conversion During Assignment**

	The above lines illustrate a feature of OpenPNM that is worth pointing out now.  All pores need to have a diffusivity value associated with them; however, we often want to assign the same value to every pore.  If you assign a scalar value to any property in OpenPNM it will automatically be converted to a vector of the appropriate length (either *Np* or *Nt* long).  This is explained in more detail :ref:`here<inner_workings>`.

The code block below illustrate how to define a *Phase* object to represent Air using 'pore-scale models'. Some of the models require various input parameters.  For instance, consider the Fuller model, which requires the molecular mass and diffusion volume of the species in the mixture.  More importantly, the Fuller model also includes temperature, meaning that if temperature of the phase changes, then the model can be re-run to regenerate the diffusivity at the new temperature.  The Fuller model code assumes that the temperature for the *Phase* can be found in ``'pore.temperature'``.  It's possible to customize these default property names as outlined :ref:`here<customizing>`.  To use the available thermophysical property models that are included with OpenPNM, import the *Phase* models library.

>>> from OpenPNM.Phases import models as fm
>>> air.add_model(propname='pore.diffusivity',
...                model=fm.diffusivity.fuller,
...                MA=0.03199,
...                MB=0.0291,
...                vA=16.3,
...                vB=19.7)
>>> air.add_model(propname='pore.viscosity',
...               model=fm.viscosity.reynolds,
...               uo=0.001,
...               b=0.1)
>>> air.add_model(propname='pore.molar_density',
...               model=fm.molar_density.ideal_gas,
...               R=8.314)

===============================================================================
Create Pore Scale Physics Objects
===============================================================================

We are still not ready to perform any simulations.  The last step is to define the desired pore scale physics, which defines how the phase and geometrical properties interact.  A classic example of this is the Washburn equation which predicts the capillary pressure required to push a non-wetting fluid through a capillary of known size.  Because the *Physics* object defines the interaction of a *Phase* with the *Geometry*, it is necessary to build one *Physics* object for each intersection between *Geometry* and *Phase* objects:

>>> phys_water = OpenPNM.Physics.GenericPhysics(network=pn,
...                                             phase=water,
...                                             geometry=geom)
>>> phys_air = OpenPNM.Physics.GenericPhysics(network=pn,
...                                           phase=air,
...                                           geometry=geom)

*Physics* objects do not require the specification of which ``pores`` and ``throats`` where they apply.  This assignment is implied by the passing of a ``geometry`` object, which has already been assigned to specific locations.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Add Desired Methods to Physics Objects
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

As with *Phase* and *Geometry* objects, the next steps are to load the model library from the *Physics* submodule and then to build-up the bare objects with the desired models:

>>> from OpenPNM.Physics import models as pm
>>> phys_water.add_model(propname='throat.capillary_pressure',
...                      model=pm.capillary_pressure.washburn)
>>> phys_water.add_model(propname='throat.hydraulic_conductance',
...                      model=pm.hydraulic_conductance.hagen_poiseuille)
>>> phys_air.add_model(propname='throat.diffusive_conductance',
...                    model=pm.diffusive_conductance.bulk_diffusion)
>>> phys_air.add_model(propname='throat.hydraulic_conductance',
...                    model=pm.hydraulic_conductance.hagen_poiseuille)

===============================================================================
Run some simulations
===============================================================================

Finally, it is now possible to run some simulations.  The code below estimates the effective diffusivity through the network by applying a concentration gradient across and calculating the flux.  This starts by creating a FickianDiffusion *Algorithm*, which is pre-defined in OpenPNM:

>>> alg = OpenPNM.Algorithms.FickianDiffusion(network=pn,phase=air)

Next the boundary conditions are applied using the ``set_boundary_conditions`` method.  In this case the boundary conditions are applied to the ``'left'`` and ``'right'`` of the cubic domain.

>>> BC1_pores = pn.pores('right')
>>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
>>> BC2_pores = pn.pores('left')
>>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.4, pores=BC2_pores)

.. note:: **Pore and Throat Labels**

	Note how the ``pores`` method was used to extract pore numbers based on the labels ``'left'`` and ``'right'``.  It's possible to add your own labels to the simulations to allow quick access to special sets of pores.  This is outlined :ref:`here<inner_workings>`.

To actually run the algorithm use the ``run`` method.  This builds the coefficient matrix from the existing values of diffusive conductance, and inverts the matrix to solve for concentration in each pores.

>>> alg.run()
>>> alg.return_results()
>>> Deff = alg.calc_eff_diffusivity()

===============================================================================
Visualise the Results
===============================================================================
We can now visualise our network and simulation results.  OpenPNM does not support native visualization, so data must be exported to a file for exploration in another program such as any of the several VTK front ends (i.e. Paraview).

.. code-block:: python

	ctrl.export(pn)

This creates a *net.vtp* file in the active directory, which can be loaded from ParaView. For a quick tutorial on the use of Paraview with OpenPNM data, see :ref:`Using Paraview<paraview_example>`.

To save an incomplete simulation for later work, the **Controller** object can be used to save the entire workspace (i.e. all simulations) using ``ctrl.save()``, or just the simulation of interest using ``ctrl.save_simulation(pn)``.
