.. _getting_started:

###############################################################################
Tutorial 1 of 3: Getting Started with OpenPNM
###############################################################################

This tutorial is intended to show the basic outline of how OpenPNM works, and skips many of the useful and powerful features of the package.  So if you find yourself asking "why is this step so labor intensive" it's probably because this tutorial deliberately simplifies some steps to provide a more gentle introduction.  The second and third tutorials of this User-Guide dive into the code more deeply, but those features are best appreciated once the basics are understood.

**Learning Objectives**

1. A brief introduction to main OpenPNM objects
2. Generate a standard cubic **Network** topology
3. Learn some handy tools for working with networks
4. Calculate geometrical properties and assign to a **Geometry** object
5. Calculate thermophysical properties and assign to a **Phase** object
6. Define pore-scale physics, and assign transport parameters to a **Physics** object
7. Run a permeability simulation using the pre-defined **Algorithm**

===============================================================================
The Core Objects
===============================================================================

OpenPNM employs 5 main objects which each store and manage a different type of data.  These are listed below:

=============  ====
Name           Role
=============  ====
**Network**    Manages the topological data such as pore locations and pore-to-pore connections.
**Geometry**   Manages the geometrical properties such as pore diameter and throat length.
**Phase**      Manages the thermophysical values such as temperature and viscosity
**Physics**    Manage the pore-scale transport parameters such as hydraulic conductance
**Algorithm**  Contain the algorithms that use the data stored on the other objects, such as diffusion or drainage
=============  ====

Each of the above objects is a *subclass* of the Python *dictionary* or *dict*, which is a very general storage container that allows values to be accessed by a name:

.. code-block:: python

	>>> foo = dict()
	>>> foo['bar'] = 1
	>>> foo['baz'] = [1, 2, 3]
	>>> sorted(foo.keys())

The *dict* does not offer much functionality aside from this rather basic storage behavior.  OpenPNM extends the *dict* to have functionality specifically suited for dealing with pore network data.  Awareness of these 5 main object types and an idea how the *dict* functions is sufficient for this tutorial.  The interested reader is referred to the detailed reference on the :ref:`object_relationship` for more information.

===============================================================================
Build a Cubic Network
===============================================================================

Start by importing OpenPNM and the Scipy package:

.. code-block:: python

	>>> import OpenPNM
	>>> import scipy as sp

Next, generate a **Network** by choosing the desired topology (e.g. cubic), then create an *instance* with the desired parameters:

.. code-block:: python

	>>> pn = OpenPNM.Network.Cubic(shape=[4, 3, 1], spacing=0.0001)

This generates a topological network using the *Cubic* class, and stores it in variable ``pn``.  This network contains pores at the correct spatial positions and connections between the pores according the cubic topology.  The ``shape`` argument specifies the number of pores in the [X, Y, Z] directions of the cube.  Networks in OpenPNM are always 3D dimensional, meaning that a 2D or "flat" network is still 1 layer of pores "thick" so [X, Y, Z] = [20, 10, 1], thus ``pn`` in this tutorial is 2D which is easier for visualization.  The ``spacing`` argument controls the center-to-center distance between pores and it can be a scalar or vector (i.e. [0.0001, 0.0002, 0.0003]).  `Using Paraview for Visualization`__ gives the network shown below:

.. image:: http://i.imgur.com/ScdydO9.png

.. note:: **Units in OpenPNM, or lack thereof**

	Although OpenPNM does not currently have a dimensional units system, it is *strongly* recommend using SI throughout.  We have investigated several packages but none are fully satisfactory...yet.

-------------------------------------------------------------------------------
Handy Tools for Working with Networks
-------------------------------------------------------------------------------

**Network** objects have numerous methods or functions that can be used to query their topological properties:

.. code-block:: python

	>>> pn.find_neighbor_pores(pores=[1])  # Find neighbors of pore 1
	array([0, 2, 4])
	>>> pn.find_neighbor_throats(pores=[1, 2])  # Find throats connected to pores 1 and 2
	array([ 0,  1,  9, 10])

There are several more such topological query methods available on **Network** objects such as ``find_nearby_pores``, ``find_connecting_throat`` and ``find_clusters``.  For more information these tools see the :ref:`topology`.

===============================================================================
Initialize and Build a Geometry Object
===============================================================================

The **Network** ``pn`` does not contain any information about pore and throat sizes at this point.  The next step, then, is to create a **Geometry** object to manage the geometrical properties.

.. code-block:: python

	>>> geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=pn.Ps,
	...                                         throats=pn.Ts)

This statement contains three arguments: ``network`` tells the **Geometry** object which **Network** it is associated with, and  ``pores`` and ``throats`` indicate the locations in the **Network** where this **Geometry** object will apply.  In this simple tutorial it is *all* pores and throats, but there are many cases where different regions of the network  have different geometrical properties, so multiple **Geometry** objects can be created, but this is a subject for the `intermediate tutorial <intermediate_usage>`_.

-------------------------------------------------------------------------------
Add Desired Size Information
-------------------------------------------------------------------------------

This freshly instantiated **Geometry** object ``geom`` contains no geometric properties as yet.  For this tutorial we'll use the direct assignment of static values (See the `intermediate tutorial <intermediate_usage>`_ for more details).

Let's start by assigning diameters to each pore from a random distribution, spanning 0 um to 100 um.  The upper limit arises because the ``spacing`` of the **Network** was set to 100 [um], so pore diameters exceeding 100 um might overlap with their neighbors.

.. code-block:: python

	>>> geom['pore.diameter'] = sp.rand(pn.Np)*0.0001

This creates an array of random numbers between and 0.0001 that is *Np-long*, meaning each pore is assigned a unique random number.

For throat diameters, we want them to always be smaller than the two pores which it connects to maintain physical consistency. This requires understanding a little bit about how OpenPNM stores network topology.  `A detailed explanation is given elsewhere <topology>`_.  Consider the following:

.. code-block:: python

	>>> P12 = pn['throat.conns']  # An Nt x 2 list of pores on the end of each throat
	>>> D12 = geom['pore.diameter'][P12]  # An Nt x 2 list of pore diameters
	>>> Dt = sp.amin(D12, axis=1)  # An Nt x 1 list of the smaller pore from each pair
	>>> geom['throat.diameter'] = Dt

Let's dissect the above lines.  Firstly, ``P12`` is a direct copy of the **Network's** ``'throat.conns'`` array, which contains the indices of the pore-pair connected by each throat.  Next, this *Nt-by-2* array is used to index into the ``'pore.diameter'`` array, resulting in another *Nt-by-2* array containing the diameters of the pores on each end of a throat.  Finally, the Scipy function ``amin`` is used to find the minimum diameter of each pore-pair by specifying the ``axis`` argument as 1, and the resulting *Nt-by-1* array is assigned to ``geom['throat.diameter']``.  This trick of using ``'throat.conns'`` to index into a pore property array is commonly used in OpenPNM and you should have a second look at the above code to understand it fully.

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

The simulation is now topologically and geometrically defined.  It has pore coordinates, pore and throat sizes and so on.  In order to perform any simulations it is necessary to define **Phase** objects that represent the fluids in the simulations:

.. code-block:: python

	>>> air = OpenPNM.Phases.GenericPhase(network=pn)
	>>> water = OpenPNM.Phases.GenericPhase(network=pn)

``pn`` is passed as an argument because **Phases** must know to which **Network** they belong.  Also, note that ``pores`` and ``throats`` are NOT specified; this is because **Phases** are mobile and can exist anywhere or everywhere in the domain, so providing specific locations does not make sense.  Algorithms for dynamically determining actual phase distributions are discussed later.

-------------------------------------------------------------------------------
Add Desired Thermophysical Properties
-------------------------------------------------------------------------------

Now it is necessary to fill these two **Phase** objects with the desired thermophysical properties.  The most basic means is to simply assign static values as follows:

.. code-block:: python

		>>> water['pore.temperature'] = 298.0
		>>> water['pore.viscosity'] = 0.001
		>>> air['pore.temperature'] = 298.0
		>>> air['pore.viscosity'] = 0.0000173

OpenPNM includes a framework for calculating these type of properties from models and correlations, but this is covered in the `intermediate tutorial <intermediate_usage>`_.

    | **Scalar to Vector Conversion During Assignment**: The above lines illustrate a feature of OpenPNM that is worth pointing out now.  All pores need to have a diffusivity value associated with them; however, we often want to assign the same value to every pore.  If you assign a scalar value to any property in OpenPNM it will automatically be converted to a vector of the appropriate length (either *Np* or *Nt* long).

===============================================================================
Create Physics Objects
===============================================================================

We are still not ready to perform any simulations.  The last step is to define the desired pore scale physics models, which dictate how the phase and geometrical properties interact.  A classic example of this is the Hagen-Poiseuille equation for fluid flow through a throat to predict the flow rate as a function of the pressure drop.  The flow rate is proportional to the geometrical size of the throat (radius and length) as well as properties of the fluid (viscosity).  It follows that this calculation needs to be performed once for each phase of interest since each has a different viscosity.  This is accomplished by define a **Physics** object for each *Phase*:

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

The viscosity of the **Phases** was only defined in the pores; however, the hydraulic conductance must be calculated for each throat.  There are several options: (1) use a hard-coded scalar value in the calculation, (2) assign ``'throat.viscosity'`` to each phase or (3) use interpolation to estimate throat viscosity as an average of the values in the neighboring pores.  The third option is suitable when there is a distribution of temperatures throughout the network and therefore viscosity changes as well, and OpenPNM provides tools for this which are discussed later.  In the present case as simple scalar value is sufficient:

.. code-block:: python

	>>> mu_w = 0.001
	>>> phys_water['throat.hydraulic_conductance'] = 3.14159*R**4/(8*mu_w*L)
	>>> mu_a = 0.0000173
	>>> phys_air['throat.hydraulic_conductance'] = 3.14159*R**4/(8*mu_a*L)

Note that both of these calculations use the same geometrical properties (``R`` and ``L``) but different phase properties (``mu_w`` and ``mu_a``).  This is why a new **Physics** object is required for each **Phase** that is added.

===============================================================================
Create an Algorithm Object for Performing a Permeability Simulation
===============================================================================

Finally, it is now possible to run some simulations.  The code below estimates the permeability through the network by applying a pressure gradient across and calculating the flux.  This starts by creating a **StokesFlow** algorithm, which is pre-defined in OpenPNM:

.. code-block:: python

	>>> alg = OpenPNM.Algorithms.StokesFlow(network=pn, phase=air)

Like all the above objects, **Algorithms** must be assigned to a **Network** via the ``network`` argument.  This algorithm is also associated with a **Phase** object, in this case ``air``, which dictates which pore-scale **Physics** properties to use (recall that ``phys_air`` was associated with ``air``).

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

The results (``'pore.pressure'``) are held within the ``alg`` object and must be explicitly returned to the ``air`` object by the user if they wish to use these values in a subsequent calculation.  The point of this data containment is to prevent unwanted overwriting of data.  Each algorithm has a method called ``return_results`` which places the pertinent values back onto the appropriate **Phase** object.

.. code-block:: python

	>>> alg.return_results()

Using Paraview for Visualization, the resulting pressure gradient across the network can be seen:

.. image:: http://i.imgur.com/8aVaH1S.png

.. _Github_Paraview_Tutorial:  https://github.com/PMEAL/OpenPNM-Examples/blob/master/IO_and_Visualization/paraview.md

__ Github_Paraview_Tutorial_
