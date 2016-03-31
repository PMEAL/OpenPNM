.. _getting_started:

###############################################################################
Tutorial 1 of 3: Getting Started with OpenPNM
###############################################################################

This tutorial is intended to show the basic outline of how OpenPNM works, and necessarily skips many of the useful and powerful features of the package.  So if you find yourself asking "why is this step so labor intensive" it's probably because this tutorial deliberately simplifies some steps to provide a more gentle introduction.  The second and third tutorials of this User-Guide dive into the code more deeply, but those features are best appreciated once the basics are understood.

.. note:: Learn Python the Hard Way

	If you're looking for an excellent online, interactive tutorial on the general use of Python, the `Learn Python the Hard Way <http://learnpythonthehardway.org/book/>`_ book and website are excellent.


**Learning Objectives**

#. Grasp the main OpenPNM objects and their roles
#. Generate a standard cubic **Network** topology
#. Learn some handy tools for working with objects and networks in particular
#. Calculate geometrical properties and assign to a **Geometry** object
#. Perform the typical data read and write operations
#. Use the network topology storage scheme to perform some calculations
#. Calculate thermophysical properties and assign to a **Phase** object
#. Define pore-scale physics, and assign transport parameters to a **Physics** object
#. Run a permeability simulation using the pre-defined **Algorithm**

===============================================================================
The Main OpenPNM Objects
===============================================================================

OpenPNM employs 5 main objects which each store and manage a different type of data.  These are listed below:

=============  ====
Name           Role
=============  ====
**Network**    Manages the topological data such as pore locations and pore-to-pore connections
**Geometry**   Manages the geometrical properties such as pore diameter and throat length
**Phase**      Manages the thermophysical properties such as temperature and viscosity
**Physics**    Manages the pore-scale transport parameters such as hydraulic conductance
**Algorithm**  Contains the algorithms that use the data from other objects to perform simulations, such as diffusion or drainage
=============  ====

Each of the above objects is a *subclass* of the Python *dictionary* or *dict*, which is a very general storage container that allows values to be accessed by a name:

.. code-block:: python

	>>> foo = dict()  # Create empty dict
	>>> foo['bar'] = 1  # Add new entry called 'bar'
	>>> foo['baz'] = [1, 2, 3]
	>>> foo['bar']  # Retrieve entries by name
	1
	>>> sorted(foo.keys())  # Inspect all entries
	['bar', 'baz']

A more detailed tutorial on dictionaries `can be found here <http://learnpythonthehardway.org/book/ex39.html>`_.  The *dict* does not offer much functionality aside from this basic storage, and it is in fact meant to be extended.  OpenPNM extends the *dict* to have functionality specifically suited for dealing with OpenPNM data.  Awareness of these 5 main object types and a familiarity with the *dict* syntax is sufficient for this tutorial, but a more information can be found in the :ref:`class_hierarchy`.

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

This generates a topological network using the *Cubic* class, and stores it in variable ``pn``.  This network contains pores at the correct spatial positions and connections between the pores according the cubic topology.  The ``shape`` argument specifies the number of pores in the [X, Y, Z] directions of the cube.  Networks in OpenPNM are always 3D dimensional, meaning that a 2D or "flat" network is still 1 layer of pores "thick" so [X, Y, Z] = [20, 10, 1], thus ``pn`` in this tutorial is 2D which is easier for visualization.  The ``spacing`` argument controls the center-to-center distance between pores and it can be a scalar or vector (i.e. [0.0001, 0.0002, 0.0003]).

-------------------------------------------------------------------------------
Tools for Inspecting Object Properties
-------------------------------------------------------------------------------
As mentioned above, each of the main objects in OpenPNM are Python *dicts* with a variety of methods and functions added that work specifically on OpenPNM data.  All of the main objects have methods for querying basic properties, like the number of pores or throats, which properties have been defined, and so on.

.. code-block:: python

	>>> pn.num_pores()
	12
	>>> pn.Np  # Shortcut to get number of pores
	12
	>>> pn.num_throats()
	17
	>>> pn.Nt
	17
	>>> pn.props()
	['pore.coords', 'pore.index', 'throat.conns']

More information about these various functions is given in :ref:`class_hierarchy`.

-------------------------------------------------------------------------------
Tools for Querying Network Topology
-------------------------------------------------------------------------------

In addition to the general methods for inspecting properties mentioned above, **Network** objects have additional functionality for performing queries on their topological data:

.. code-block:: python

	>>> pn.find_neighbor_pores(pores=[1])  # Find neighbors of pore 1
	array([0, 2, 4])
	>>> pn.find_neighbor_throats(pores=[1, 2])  # Find throats connected to pores 1 and 2
	array([ 0,  1,  9, 10])

There are several more such topological query methods available on **Network** objects such as ``find_nearby_pores``, ``find_connecting_throat`` and ``find_clusters``.  For more information on these tools see the :ref:`topology`.

-------------------------------------------------------------------------------
Exporting Data for Visualization
-------------------------------------------------------------------------------
OpenPNM does not offer it's own visualization tools, as there are already many excellent options available.  The workflow for visualization is to output the simulation data to a standard file format for use in a program like `Paraview <http://www.paraview.org>`_.  The most convenient way to export data is to use the ``export_data`` method in the main OpenPNM namespace:

.. code-block:: python

	>>> OpenPNM.export_data(network=pn, filename='test', fileformat='VTK')

This creates a file called *test.vtp* in the current working directory.  Note that *VTK* stands for Visualization Toolkit, and is the general name for this type of file, but the file has a *vtp* extension which is a specific type of *VTK* file.  Opening this file in Paraview gives the following result:

.. image:: http://i.imgur.com/ScdydO9.png

For help using Paraview, see the `Example in the OpenPNM-Example collection <https://github.com/PMEAL/OpenPNM-Examples/blob/master/IO_and_Visualization/paraview.md>`_

===============================================================================
Initialize and Build a Geometry Object
===============================================================================

The **Network** ``pn`` does not contain any information about pore and throat sizes at this point.  The next step, then, is to create a **Geometry** object to manage the geometrical properties.

.. code-block:: python

	>>> geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=pn.Ps,
	...                                         throats=pn.Ts)

This statement contains three arguments: ``network`` tells the **Geometry** object which **Network** it is associated with, and  ``pores`` and ``throats`` indicate the locations in the **Network** where this **Geometry** object will apply.  In this  tutorial ``geom`` applies to *all* pores and throats, but there are many cases where different regions of the network have different geometrical properties, so OpenPNM allows multiple **Geometry** objects to be created for managing the data in each region, but this is a subject for :ref:`intermediate_usage`.

-------------------------------------------------------------------------------
Add Desired Size Information
-------------------------------------------------------------------------------

This freshly instantiated **Geometry** object (``geom``) contains no geometric properties as yet.  For this tutorial we'll use the direct assignment of manually calculated values.

Let's start by assigning diameters to each pore from a random distribution, spanning 0 um to 100 um.  The upper limit matches the ``spacing`` of the **Network** was set to 100 [um], so pore diameters exceeding 100 um might overlap with their neighbors.

.. code-block:: python

	>>> geom['pore.diameter'] = sp.rand(pn.Np)*0.0001

This creates an array of random numbers between 0 and 0.0001 that is *Np*-long, meaning each pore is assigned a unique random number.

The previous line illustrates a key point about data storage rules in OpenPNM.  Note that the array name started with ``'pore.'``.  All dictionary entries must start with either ``'pore.'`` or ``'throat.'``.  The reason for this is that OpenPNM forces arrays to be of the appropriate length (either *Nt* or *Np* long), which it infers from the name of the array.  Attempts to write arrays of the wrong length are blocked:

.. code-block:: python

	>>> geom['foo'] = sp.ones(pn.Np)  # Will result in an exception
	>>> geom['pore.foo'] = sp.ones(pn.Np - 2)  # Will result in an error message
	>>> geom['throat.foo'] = sp.one(pn.Np)  # Also gives an error message

Returning to the definition of **Geometry** properties, we want the throat diameters to always be smaller than the two pores which it connects to maintain physical consistency. This requires understanding a little bit about how OpenPNM stores network topology.  Consider the following:

.. code-block:: python

	>>> P12 = pn['throat.conns']  # An Nt x 2 list of pores on the end of each throat
	>>> D12 = geom['pore.diameter'][P12]  # An Nt x 2 list of pore diameters
	>>> Dt = sp.amin(D12, axis=1)  # An Nt x 1 list of the smaller pore from each pair
	>>> geom['throat.diameter'] = Dt

Let's dissect the above lines.  Firstly, ``P12`` is a direct copy of the **Network's** ``'throat.conns'`` array, which contains the indices of the pore-pair connected by each throat.  Next, this *Nt-by-2* array is used to index into the ``'pore.diameter'`` array, resulting in another *Nt-by-2* array containing the diameters of the pores on each end of a throat.  Finally, the Scipy function ``amin`` is used to find the minimum diameter of each pore-pair by specifying the ``axis`` argument as 1, and the resulting *Nt-by-1* array is assigned to ``geom['throat.diameter']``.  This trick of using ``'throat.conns'`` to index into a pore property array is commonly used in OpenPNM and you should have a second look at the above code to understand it fully.  Refer to :ref:`topology` for a full discussion.

We must still specify the remaining geometrical properties of the pores and throats. Since we're creating a "Stick-and-Ball" geometry, the sizes are calculated from the geometrical equations for spheres and cylinders.

For pore volumes, assume a sphere:

.. code-block:: python

	>>> Rp = geom['pore.diameter']/2
	>>> geom['pore.volume'] = (4/3)*3.14159*(Rp)**3

The length of each throat is the center-to-center distance between pores, minus the radius of each of two neighboring pores.

.. code-block:: python

	>>> C2C = 0.0001  # The center-to-center distance between pores
	>>> Rp12 = Rp[pn['throat.conns']]
	>>> geom['throat.length'] = C2C - sp.sum(Rp12, axis=1)

The volume of each throat is found assuming a cylinder:

.. code-block:: python

    >>> Rt = geom['throat.diameter']/2
    >>> Lt = geom['throat.length']
    >>> geom['throat.volume'] = 3.14159*(Rt)**2*Lt

The basic geometrical properties of the network are now defined.  The **Geometry** class possess a method called ``plot_histograms`` that produces a plot of the most pertinent geometrical properties.  The following figure doesn't look very good since the network in this example has only 12 pores, but the utility of the plot for quick inspection is apparent.

.. image:: http://i.imgur.com/xkK1TYf.png

===============================================================================
Create Phases
===============================================================================

The simulation is now topologically and geometrically defined.  It has pore coordinates, pore and throat sizes and so on.  In order to perform any simulations it is necessary to define **Phase** objects that represent the fluids in the simulations:

.. code-block:: python

	>>> water = OpenPNM.Phases.GenericPhase(network=pn)

``pn`` is passed as an argument because **Phases** must know to which **Network** they belong.  Also, note that ``pores`` and ``throats`` are NOT specified; this is because **Phases** are mobile and can exist anywhere or everywhere in the domain, so providing specific locations does not make sense.  Algorithms for dynamically determining actual phase distributions are discussed later.

-------------------------------------------------------------------------------
Add Desired Thermophysical Properties
-------------------------------------------------------------------------------

Now it is necessary to fill these two **Phase** objects with the desired thermophysical properties.  The most basic means is to simply assign static values as follows:

.. code-block:: python

		>>> water['pore.temperature'] = 298.0
		>>> water['pore.viscosity'] = 0.001

The above code block highlight another key feature of data storage in OpenPNM.  When a scalar value is written to an object it is extended to a vector of the appropriate length (either *Np* or *Nt*) depending on the name of the array.  Although this is slightly wasteful of memory, it vastly simplifies data access since all values are explicitly defined on every pore and throat:

.. code-block:: python

	>>> water.Np
	12
	>>> len(water['pore.temperature'])
	12
	>>> water['pore.temperature'][10]
	298.0

OpenPNM includes a framework for calculating these type of properties from models and correlations, but this is covered in :ref:`intermediate_usage`.

===============================================================================
Create Physics Objects
===============================================================================

We are still not ready to perform any simulations.  The last step is to define the desired pore-scale physics models, which dictate how the phase and geometrical properties interact.  A classic example of this is the Hagen-Poiseuille equation for fluid flow through a throat to predict the flow rate as a function of the pressure drop.  The flow rate is proportional to the geometrical size of the throat (radius and length) as well as properties of the fluid (viscosity):

.. code-block:: python

	>>> phys_water = OpenPNM.Physics.GenericPhysics(network=pn, phase=water,
	...                                             geometry=geom)

**Physics** objects do not require the specification of which ``pores`` and ``throats`` where they apply, since this information is implied by the ``geometry`` argument which was already assigned to specific locations.

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

===============================================================================
Create an Algorithm Object for Performing a Permeability Simulation
===============================================================================

Finally, it is now possible to run some simulations.  The code below estimates the permeability through the network by applying a pressure gradient across and calculating the flux.  This starts by creating a **StokesFlow** algorithm, which is pre-defined in OpenPNM:

.. code-block:: python

	>>> alg = OpenPNM.Algorithms.StokesFlow(network=pn, phase=water)

Like all the above objects, **Algorithms** must be assigned to a **Network** via the ``network`` argument.  This algorithm is also associated with a **Phase** object, in this case ``water``, which dictates which pore-scale **Physics** properties to use (recall that ``phys_water`` was associated with ``water``).

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

To determine the permeability coefficient, we must invoke Darcy's law: Q = KA/uL(Pin - Pout).  Everything in this equation is known except for the volumetric flow rate Q.  The **StokesFlow** algorithm possesses a ``rate`` method that calculates the rate of a quantity leaving a specified set of pores:

.. code-block:: python

	>>> Q = alg.rate(pores='top')
	>>> A = 0.0001*3*1  # Cross-sectional area for flow
	>>> L = 0.0001*4  # Length of flow path
	>>> del_P = 101325  # Specified pressure gradient
	>>> K = Q*mu_w*L/(A*del_P)

The results (``'pore.pressure'``) are held within the ``alg`` object and must be explicitly returned to the ``air`` object by the user if they wish to use these values in a subsequent calculation.  The point of this data containment is to prevent unintentional overwriting of data.  Each algorithm has a method called ``return_results`` which places the pertinent values back onto the appropriate **Phase** object.

.. code-block:: python

	>>> alg.return_results()

Using Paraview for Visualization, the resulting pressure gradient across the network can be seen:

.. image:: http://i.imgur.com/8aVaH1S.png
