.. _intermediate_usage:

###############################################################################
Tutorial 2: Digging Deeper with OpenPNM
###############################################################################

This tutorial will follow the same outline as the :ref:`getting started tutorial <getting_started>`, but will dig a little bit deeper at each step to reveal the more advanced features and usages of OpenPNM.  Be sure you've done and understood that tutorial before attempting this one.

As usual, start by importing the OpenPNM package and the Scipy package which is always handy:

>>> import OpenPNM
>>> import scipy as sp

===============================================================================
Building a Cubic Network
===============================================================================

Let's generate a cubic network but with a different connectivity:

>>> pn = OpenPNM.Network.Cubic(shape=[20, 20, 5],
...                            spacing=0.0001,
...                            connectivity=8)

This **Network** has pores distributed in a cubic lattice, but connected to diagonal neighbors due to the ``connectivity`` being set to 8 (the default is 6).  The various options are outlined in the *Cubic* class's documentation which can be viewed with the object inspector.

OpenPNM includes several other classes for generating networks including random topology based on Delaunay tessellations (**Delaunay**), and a class for importing networks from external code such as network extractions (**Import**).

===============================================================================
Initialize and Build a Geometry Object
===============================================================================

In this tutorial we will make materials that has different geometrical properties in different regions.  This will demonstrate the motivation behind separating the **Geometry** properties from the **Network** topology.  Let's say that the pores on the top and bottom surfaces are substantially smaller than the internal pores.  We need to create one **Geometry** object to manage the top and bottom pores, and a second to manage the remaining internal pores:

>>> Ps = pn.pores(['top', 'bottom'])
>>> geom1 = OpenPNM.Geometry.GenericGeometry(network=pn, pores=Ps, name='surface')
>>> Ps = pn.pores(['top', 'bottom'], mode='not')
>>> geom2 = OpenPNM.Geometry.GenericGeometry(network=pn, pores=Ps, name='core')

The above statements result in two distinct **Geometry** objects each applying to different regions of the full network domain.  As we shall see, this simplifies the data management in some important ways.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Add Desired Properties to Geometry
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This freshly instantiated *Geometry* object contains no geometric properties as yet because we chose to use the *GenericGeometry* class.  There are several other 'pre-written' classes available in *Geometry* module, but for this example we will explicity make a 'stick-and-ball' geometry, which as the name suggests treats the pores as spheres and the throats as cylinders.

-------------------------------------------------------------------------------
Direct Assignment of Static Values
-------------------------------------------------------------------------------

Let's start by assiging diameters to each pore from a random distribution, spanning 10 um to 100 um.  The upper limit is arise because the ``spacing`` of the *Network* was set to 100 [um], so pore diameters exceeding 100 um might overlap with neighbors.  The lower limit is merely to avoid vanishingly small pores.

>>> import scipy as sp  # Import the Scipy package
>>> geom['pore.diameter'] = 0.00001 + sp.rand(pn.Np)*0.00099

This creates a ND-array of random numbers (between 0.00001 and 0.0001) that is *Np* long, meaning each pore is assigned a unique random number.

For throat diameter, we want them to always be smaller than the two pores which it connects to maintain physical consistency. This requires explaining how OpenPNM stores network topology.

>>> P12 = pn['throat.conns']  # An Nt x 2 list of pores on the end of each throat
>>> D12 = geom['pore.diameter'][P12]  # An Nt x 2 list of pore diameters
>>> Dt = sp.amin(D12, axis=1)  # An Nt x 1 list of the smaller pore from each pair
>>> geom['throat.diameter'] = Dt

Let's disect the above lines.  Firstly, P12 is a direct copy of the Network's \'throat.conns\' array, which contains the indices of the pore pair connected by each throat.  Next, this *Nt-by-2* array is used to index into the \'pore.diameter'\ array, resulting in another *Nt-by-2* array containing the diameters of the pores connected by each throat.  Finally, the Scipy function ``amin`` is used to find the minimum diameter of each pore pair by specifying the ``axis`` keyword as 1, and the resulting *Nt-by-1* array is assigned to ``geom['throat.diameter']``.

Finally, we must specify the remaining geometrical properties of the pores and throats. Since we're creating a 'stick-and-ball' geometry, the sizes are calculated from the geometrical equations for spheres and cylinders.

For pore volumes, assume a sphere:

>>> Rp = geom['pore.diameter']/2
>>> geom['pore.volume'] = (4/3)*3.14159*(Rp)**3

The length of each throat is the center-to-center distance between pores, minus the radius of each of two neighbor pores.

>>> C2C = 0.0001  # The center-to-center distance between pores
>>> Rp12 = Rp[pn['throat.conns']]
>>> geom['throat.length'] = C2C - sp.sum(Rp12, axis=1)

The volume of each throat is found assuming a cylinder:

>>> Rt = geom['throat.diameter']/2
>>> Lt = geom['throat.length']
>>> geom['throat.volume'] = 3.14159*(Rt)**2*Lt

The basic geometrical properties of the network are now defined.

===============================================================================
Create Phases
===============================================================================

The simulation is now topologically and geometrically complete.  It has pore coordinates, pore and throat sizes and so on.  In order to perform any simulations it is necessary to define *Phase* objects that represent the fluids in the simulations:

>>> air = OpenPNM.Phases.GenericPhase(network=pn, name='air')
>>> water = OpenPNM.Phases.GenericPhase(network=pn, name='water')

``pn`` is passed as an argument because *Phases* must know to which *Network* they belong.  Also, note that ``pores`` and ``throats`` are NOT specified; this is because *Phases* are mobile and can exist anywhere or everywhere in the domain, so providing specific locations does not make sense.  Algorithms for dynamically determining actual phase distributions are discussed later.

.. note:: **Naming Objects**

	The above two lines also include a ``name`` argument.  All objects in OpenPNM can be named in this way if desired, however, if no name is given one will be generated.  The point of the name is to allow easy identification of an object at the command line, using the ``name`` attribute (``air.name``).  Objects can be renamed, so if you wish to override a default name simply use ``air.name = 'air'``.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Add Desired Properties to Phases
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Now it is necessary to fill these two *Phase* objects with the desired thermophysical properties.  The most basic means is to simply assign static values as follows:

>>> water['pore.temperature'] = 298.0
>>> water['pore.viscosity'] = 0.001
>>> air['pore.temperature'] = 298.0
>>> air['pore.viscosity'] = 0.0000173

OpenPNM includes a framework for calculating these type of properties from models and correlations, but this is beyond the aim of the present introductory tutorial.

.. note:: **Scalar to Vector Conversion During Assignment**

	The above lines illustrate a feature of OpenPNM that is worth pointing out now.  All pores need to have a diffusivity value associated with them; however, we often want to assign the same value to every pore.  If you assign a scalar value to any property in OpenPNM it will automatically be converted to a vector of the appropriate length (either *Np* or *Nt* long).  This is explained in more detail :ref:`here<inner_workings>`.

===============================================================================
Create Pore Scale Physics Objects
===============================================================================

We are still not ready to perform any simulations.  The last step is to define the desired pore scale physics models, which dictates how the phase and geometrical properties interact.  A classic example of this is the Hagen-Poiseuille equation for fluid flow through a throat, which predicts the flow rate as a function of the pressure drop  The flow rate is proportional to the geometrical size of the throat (radius and length) as well as properties of the fluid (viscosity).  It follows that this calculation needs to be performed once for each phase of interest since each has a different visocity.  This is accomlished by define a *Physics* object for each *Phase*:

>>> phys_water = OpenPNM.Physics.GenericPhysics(network=pn,
...                                             phase=water,
...                                             geometry=geom)
>>> phys_air = OpenPNM.Physics.GenericPhysics(network=pn,
...                                           phase=air,
...                                           geometry=geom)

*Physics* objects do not require the specification of which ``pores`` and ``throats`` where they apply, since this information is provided by the ``geometry`` argument which has already been assigned to specific locations.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Specify Desired Pore-Scale Models
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

We need to calculate the numerical values representing our chosen pore-scale physics.  To continue with the Hagen-Poiseuille example lets calculte the hydraulic conductance of each throat in the network.  The throat radius and length are easily accessed as:

>>> R = geom['throat.diameter']/2
>>> L = geom['throat.length']

The viscosity of the *Phases* was only defined in the pores; however, the hydraulic conductance must be calculated for each throat.  There are several options: (1) use a scalar value, (2) assign \'throat.viscosity\' to each phase or (3) use interpolation to estimate throat viscosty as an average of the values in the neighboring pores.  The third option is suitable when there is a distribution of temperatures throughout the network and therefore visocity changes as well, and OpenPNM provides tools for this which are discussed later.  In the present case as simple scalar value is sufficient:

>>> mu_w = 0.001
>>> phys_water['throat.hydraulic_conductance'] = 3.14159*R**4/(8*mu_w*L)
>>> mu_a = 0.0000173
>>> phys_air['throat.hydraulic_conductance'] = 3.14159*R**4/(8*mu_a*L)

Note that both of these calcualation use the same geometrical properties (R and L) but different phase properties (mu_w and mu_a).

===============================================================================
Run Some Simulations
===============================================================================

Finally, it is now possible to run some simulations.  The code below estimates the permeabilty through the network by applying a pressure gradient across and calculating the flux.  This starts by creating a StokesFlow *Algorithm*, which is pre-defined in OpenPNM:

>>> alg = OpenPNM.Algorithms.StokesFlow(network=pn, phase=air)

Like all the above objects, algorithms must be assigned to a *Network* via the ``network`` argument.  This algorithm is also associated with a *Phase* object, in this case ``air``, which dictates which pore-scale *Physics* properties to use (recall that ``phys_air`` was associated with ``air``).

Next the boundary conditions are applied using the ``set_boundary_conditions`` method on the *Algorithm* object.  Let's apply a 1 atm pressure gradient between the left and right sides of the domain:

>>> BC1_pores = pn.pores('right')
>>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=202650, pores=BC1_pores)
>>> BC2_pores = pn.pores('left')
>>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=101325, pores=BC2_pores)

.. note:: **Pore and Throat Labels**

	Note how the ``pores`` method was used to extract pore numbers based on the labels 'left' and 'right'.  It's possible to add your own labels to simulations to allow quick access to special sets of pores.  This is outlined :ref:`here<inner_workings>`.

To actually run the algorithm use the ``run`` method.  This builds the coefficient matrix from the existing values of hydraulic conductance, and inverts the matrix to solve for pressure in each pore, and stores the results within the *Algorithm's* dictionary under \'pore.pressure'\:

>>> alg.run()

The results ('pore.pressure') are held within the ``alg`` object and must be explicitly returned to the ``air`` object by the user if they wish to use these values in a subsequent calcualation.  The point of this data containment is to prevent unwanted overwriting of data.  Each algorithm has a method called ``return_results`` which places the pertinent values back onto the appropriate *Phase* object.

>>> alg.return_results()

===============================================================================
Visualise the Results
===============================================================================
We can now visualise our network and simulation results.  OpenPNM does not support native visualization, so data must be exported to a file for exploration in another program such as any of the several VTK front ends (i.e. Paraview).

>>> OpenPNM.export(network=pn, filename='net.vtp')

This creates a *net.vtp* file in the active directory, which can be loaded from ParaView. For a quick tutorial on the use of Paraview with OpenPNM data, see :ref:`Using Paraview<paraview_example>`.
