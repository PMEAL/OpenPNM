.. _getting_started:

###############################################################################
Tutorial 1 of 3: Getting Started with OpenPNM
###############################################################################

This tutorial is intended to show the basic outline of how OpenPNM works, and necessarily skips many of the useful and powerful features of the package.  So, if you find yourself asking "why is this step so labor intensive" it's probably because this tutorial dilberately simplifies all steps to provide a more gentle introduction.  The second and third tutorials of this guide dive into the code more deeply, but those features are best appreciated once the basis are understood.  

**Learning Objectives**

1. Create a standard cubic network topology
2. Explore some handy tools methods for working with networks
3. Create a geometry object and assign geometrical properties
4. Create a phase object and assign thermo-physical properties
5. Create a physics object and calculate pore-scale transport parameters

===============================================================================
Build a Cubic Network
===============================================================================

Start by importing the OpenPNM package, and the Scipy package:

.. code-block:: python

		>>> import OpenPNM
		>>> import scipy as sp

Next, generate a **Network** by choosing the desired topology (e.g. cubic), then initialize it with the desired parameters:

.. code-block:: python

		>>> pn = OpenPNM.Network.Cubic(shape=[4, 3, 1], spacing=0.0001)

This generates a topological network and stores it in variable ``pn``.  This network contains pores at the correct spatial positions and connections between the pores according the specified topology.  The ``shape`` argument specifies the number of pores in the [X, Y, Z] directions of the cube.  Networks in OpenPNM are always 3D dimensional, meaning that a 2D or "flat" network is still 1 layer of pores "thick" so [X, Y, Z] = [20, 10, 1], so ``pn`` in this tutorial is 2D.  The ``spacing`` argument controls the center-to-center distance between pores.  Although OpenPNM does not currently have a dimensional units system, it is *strongly* recommend using SI throughout.  `Using Paraview for Visualization`_ gives the  resulting network as shown below:

.. image:: http://i.imgur.com/ScdydO9.png

-------------------------------------------------------------------------------
Handy Tools for Working with Networks
-------------------------------------------------------------------------------

The **Network** object has numerous methods that can be used to query its topological properties:

.. code-block:: python

		>>> pn.num_pores()
		12
		>>> pn.num_throats()
		17
		>>> pn.find_neighbor_pores(pores=[1])  # Find neighbors of pore 1
		array([0, 2, 4])
		>>> pn.find_neighbor_throats(pores=[1, 2])  # Find throats connected to pores 1 and 2
		array([ 0,  1,  9, 10])

There are several more such topological query method available on the object such as ``find_nearby_pores``, ``find_connecting_throat`` and ``find_clusters``.

-------------------------------------------------------------------------------
Labelling Pores and Throats
-------------------------------------------------------------------------------

Another important feature is the use of *labels* on pores and throats.  Applying a label to a set of special pores allows for easy retrieval of these pores for later use.  For instance, during the generation of a **Cubic** network, the faces are automatically labeled (i.e. 'top', 'front').  The following illustrates how to use *labels*:

.. code-block:: python

		>>> pn.labels(pores=[1])  # Find all labels applied to pore 1
		['pore.all', 'pore.bottom', 'pore.front', 'pore.internal', 'pore.top']

Note here that *pore* 1 is labelled both ``'pore.top'`` and ``'pore.bottom'`` since this is a 2D network.

The following illustrates the ``mode`` option, which controls the *set logic* to apply to the query:

.. code-block:: python

		>>> pn.pores(labels=['front', 'left'], mode='intersection')
		array([0])

Adding custom *labels* is also straight-forward.  A *label* is simply a Boolean array stored as ``'pore.some_label'`` with True values where the *label* applies. For instance, if the front-left corner of the network were of special importance and you would like to easily retrieve these pores in the future you could label them:

.. code-block:: python

		>>> Ps = pn.pores(labels=['front', 'left'], mode='intersection')
		>>> pn['pore.front_left_corner'] = False
		>>> pn['pore.front_left_corner'][Ps] = True

Note that we had to create an array for the label first filled with False values, and then assign True values in the locations where the label ``'front_left_corner'`` applies.  `Using Paraview for Visualization`_ pores labelled ``'front_left_corner'`` are colored in red:

.. image:: http://i.imgur.com/RE5DjzS.png

===============================================================================
Initialize and Build a Geometry Object
===============================================================================

The **Network** ``pn`` does not contain any information about pore and throat sizes at this point.  The next step, then, is to create a **Geometry** object to calculate the desired geometrical properties.

.. code-block:: python

		>>> geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=pn.Ps,
		...                                         throats=pn.Ts)

This statement contains three arguments: ``network`` tells the **Geometry** object which **Network** it is associated with.  ``pores`` and ``throats`` indicate the locations in the **Network** where this **Geometry** object will apply.  In this case it is all pores and throats (See the intermediate tutorial for more details).

-------------------------------------------------------------------------------
Add Desired Size Information
-------------------------------------------------------------------------------

This freshly instantiated **Geometry** object ``geom`` contains no geometric properties as yet.  For this tutorial we'll use the direct assignment of static values (See the intermediate tutorial for more details).

Let's start by assigning diameters to each pore from a random distribution, spanning 10 um to 100 um.  The upper limit arises because the ``spacing`` of the **Network** was set to 100 [um], so pore diameters exceeding 100 um might overlap with their neighbors.  The lower limit is to avoid vanishingly small pores.

.. code-block:: python

		>>> geom['pore.diameter'] = 0.00001 + sp.rand(pn.Np)*0.000099

This creates an array of random numbers (between 0.00001 and 0.0001) that is *Np-long*, meaning each pore is assigned a unique random number.

For throat diameters, we want them to always be smaller than the two pores which it connects to maintain physical consistency. This requires explaining how OpenPNM stores network topology.  Consider the following:

.. code-block:: python

		>>> P12 = pn['throat.conns']  # An Nt x 2 list of pores on the end of each throat
		>>> D12 = geom['pore.diameter'][P12]  # An Nt x 2 list of pore diameters
		>>> Dt = sp.amin(D12, axis=1)  # An Nt x 1 list of the smaller pore from each pair
		>>> geom['throat.diameter'] = Dt

Let's dissect the above lines.  Firstly, ``P12`` is a direct copy of the **Network's** ``'throat.conns'`` array, which contains the indices of the pore-pair connected by each throat.  Next, this *Nt-by-2* array is used to index into the ``'pore.diameter'`` array, resulting in another *Nt-by-2* array containing the diameters of the pores on each end of a throat.  Finally, the Scipy function ``amin`` is used to find the minimum diameter of each pore-pair by specifying the ``axis`` argument as 1, and the resulting *Nt-by-1* array is assigned to ``geom['throat.diameter']``.

We must still specify the remaining geometrical properties of the pores and throats. Since we're creating a "Stick-and-Ball" geometry, the sizes are calculated from the geometrical equations for spheres and cylinders.

For pore volumes, assume a sphere:

.. code-block:: python

		>>> Rp = geom['pore.diameter']/2
		>>> geom['pore.volume'] = (4/3)*3.14159*(Rp)**3

The length of each throat is the center-to-center distance between pores, minus the radius of each of two neighbor pores.

.. code-block:: python

		>>> C2C = 0.0001  # The center-to-center distance between pores
		>>> Rp12 = Rp[pn['throat.conns']]
		>>> geom['throat.length'] = C2C - sp.sum(Rp12, axis=1)

The volume of each throat is found assuming a cylinder:

.. code-block:: python

    >>> Rt = geom['throat.diameter']/2
    >>> Lt = geom['throat.length']
    >>> geom['throat.volume'] = 3.14159*(Rt)**2*Lt

The basic geometrical properties of the network are now defined.  The **Geometry** class possess a method called ``plot_histograms`` that produces a plot of the most pertinent geometrical properties.  The following figure doesn't look very good since our example network only has 12 pores, but the utility of the plot should be apparent.

.. image:: http://i.imgur.com/xkK1TYf.png

===============================================================================
Create Phases
===============================================================================

The simulation is now topologically and geometrically complete.  It has pore coordinates, pore and throat sizes and so on.  In order to perform any simulations it is necessary to define **Phase** objects that represent the fluids in the simulations:

.. code-block:: python

		>>> air = OpenPNM.Phases.GenericPhase(network=pn, name='air')
		>>> water = OpenPNM.Phases.GenericPhase(network=pn, name='water')

``pn`` is passed as an argument because **Phases** must know to which **Network** they belong.  Also, note that ``pores`` and ``throats`` are NOT specified; this is because **Phases** are mobile and can exist anywhere or everywhere in the domain, so providing specific locations does not make sense.  Algorithms for dynamically determining actual phase distributions are discussed later.

    | **Naming Objects**: The above two lines also include a ``name`` argument. All objects in OpenPNM can be named in this way if desired; however, if no name is given one will be generated.  The point of the name is to allow easy identification of an object at the command line, using the ``name`` attribute  (``air.name``).  Objects can be renamed, so if you wish to override a default name simply use ``air.name = 'air'``.

-------------------------------------------------------------------------------
Add Desired Thermophysical Properties
-------------------------------------------------------------------------------

Now it is necessary to fill these two **Phase** objects with the desired thermophysical properties.  The most basic means is to simply assign static values as follows:

.. code-block:: python

		>>> water['pore.temperature'] = 298.0
		>>> water['pore.viscosity'] = 0.001
		>>> air['pore.temperature'] = 298.0
		>>> air['pore.viscosity'] = 0.0000173

OpenPNM includes a framework for calculating these type of properties from models and correlations, but this is beyond the aim of the present introductory tutorial.


    | **Scalar to Vector Conversion During Assignment**: The above lines illustrate a feature of OpenPNM that is worth pointing out now.  All pores need to have a diffusivity value associated with them; however, we often want to assign the same value to every pore.  If you assign a scalar value to any property in OpenPNM it will automatically be converted to a vector of the appropriate length (either *Np* or *Nt* long).

===============================================================================
Create Physics Objects
===============================================================================

We are still not ready to perform any simulations.  The last step is to define the desired pore scale physics models, which dictates how the phase and geometrical properties interact.  A classic example of this is the Hagen-Poiseuille equation for fluid flow through a throat, which predicts the flow rate as a function of the pressure drop  The flow rate is proportional to the geometrical size of the throat (radius and length) as well as properties of the fluid (viscosity).  It follows that this calculation needs to be performed once for each phase of interest since each has a different viscocity.  This is accomplished by define a **Physics** object for each *Phase*:

.. code-block:: python

		>>> phys_water = OpenPNM.Physics.GenericPhysics(network=pn,
		...                                             phase=water,
		...                                             geometry=geom)
		>>> phys_air = OpenPNM.Physics.GenericPhysics(network=pn,
		...                                           phase=air,
		...                                           geometry=geom)

**Physics** objects do not require the specification of which ``pores`` and ``throats`` where they apply, since this information is provided by the ``geometry`` argument which has already been assigned to specific locations.

-------------------------------------------------------------------------------
Specify Desired Pore-Scale Physics Models
-------------------------------------------------------------------------------

We need to calculate the numerical values representing our chosen pore-scale physics.  To continue with the Hagen-Poiseuille example lets calculate the hydraulic conductance of each throat in the network.  The throat radius and length are easily accessed as:

.. code-block:: python

		>>> R = geom['throat.diameter']/2
		>>> L = geom['throat.length']

The viscosity of the **Phases** was only defined in the pores; however, the hydraulic conductance must be calculated for each throat.  There are several options: (1) use a scalar value, (2) assign ``'throat.viscosity'`` to each phase or (3) use interpolation to estimate throat viscosity as an average of the values in the neighboring pores.  The third option is suitable when there is a distribution of temperatures throughout the network and therefore viscosity changes as well, and OpenPNM provides tools for this which are discussed later.  In the present case as simple scalar value is sufficient:

.. code-block:: python

		>>> mu_w = 0.001
		>>> phys_water['throat.hydraulic_conductance'] = 3.14159*R**4/(8*mu_w*L)
		>>> mu_a = 0.0000173
		>>> phys_air['throat.hydraulic_conductance'] = 3.14159*R**4/(8*mu_a*L)

Note that both of these calculation use the same geometrical properties (R and L) but different phase properties (mu_w and mu_a).  This is why a new **Physics** object is required for each **Phase** that is added.

===============================================================================
Create an Algorithm Object for Performing a Permeability Simulation
===============================================================================

Finally, it is now possible to run some simulations.  The code below estimates the permeability through the network by applying a pressure gradient across and calculating the flux.  This starts by creating a StokesFlow **Algorithm**, which is pre-defined in OpenPNM:

.. code-block:: python

		>>> alg = OpenPNM.Algorithms.StokesFlow(network=pn, phase=air)

Like all the above objects, algorithms must be assigned to a **Network** via the ``network`` argument.  This algorithm is also associated with a **Phase** object, in this case ``air``, which dictates which pore-scale **Physics** properties to use (recall that ``phys_air`` was associated with ``air``).

Next the boundary conditions are applied using the ``set_boundary_conditions`` method on the **Algorithm** object.  Let's apply a 1 atm pressure gradient between the left and right sides of the domain:

.. code-block:: python

	>>> BC1_pores = pn.pores('front')
	>>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=202650,
	...                             pores=BC1_pores)
	>>> BC2_pores = pn.pores('back')
	>>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=101325,
	...                             pores=BC2_pores)

To actually run the algorithm use the ``run`` method:

.. code-block:: python

		>>> alg.run()

This builds the coefficient matrix from the existing values of hydraulic conductance, and inverts the matrix to solve for pressure in each pore, and stores the results within the **Algorithm's** dictionary under ``'pore.pressure'``.

The results ('pore.pressure') are held within the ``alg`` object and must be explicitly returned to the ``air`` object by the user if they wish to use these values in a subsequent calculation.  The point of this data containment is to prevent unwanted overwriting of data.  Each algorithm has a method called ``return_results`` which places the pertinent values back onto the appropriate **Phase** object.

.. code-block:: python

		>>> alg.return_results()

`Using Paraview for Visualization`_ , the resulting pressure gradient across the network can be seen:

.. image:: http://i.imgur.com/8aVaH1S.png

===============================================================================
Using Paraview for Visualization
===============================================================================
We can now visualize our network and simulation results.  OpenPNM does not support native visualization, so data must be exported to a file for exploration in another program such as any of the several VTK front ends (i.e. Paraview).

.. code-block:: python

		>>> OpenPNM.export_data(network=pn, filename='2D_net')

This creates a *net.vtp* file in the active directory, which can be loaded from ParaView. For a quick tutorial on the use of Paraview with OpenPNM data, see :ref:`Using Paraview<paraview_example>`.
