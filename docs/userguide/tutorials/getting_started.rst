.. _getting_started:

.. sectnum::
   :start: 1

###############################################################################
 Tutorial 1 of 3: Getting Started with OpenPNM
###############################################################################

This tutorial is intended to show the basic outline of how OpenPNM works, and necessarily skips many of the more useful and powerful features of the package.  So if you find yourself asking "why is this step so labor intensive" it's probably because this tutorial deliberately simplifies some features to provide a more gentle introduction.  The second and third tutorials of this User-Guide dive into the package more deeply, but those features are best appreciated once the basics are understood.

.. contents:: Topics Covered in this Tutorial

**Learning Objectives**

#. Introduce the main OpenPNM objects and their roles
#. Explore the way OpenPNM stores data, including network topology
#. Learn some handy tools for working with objects
#. Generate a standard cubic **Network** topology
#. Calculate geometrical properties and assign them to a **Geometry** object
#. Calculate thermophysical properties and assign to a **Phase** object
#. Define pore-scale physics and assign transport parameters to a **Physics** object
#. Run a permeability simulation using the pre-defined **Algorithm**
#. Use the package to calculate the permeability coefficient of a porous media

.. hint:: Python and Numpy Tutorials

	* OpenPNM is written in Python.  One of the best guides to learning Python is the set of Tutorials available on the `official Python website <https://docs.python.org/3.5/tutorial>`_). The web is literally overrun with excellent Python tutorials owing to the popularity and importance of the language.  The official Python website also provides `an long list of resources <https://www.python.org/about/gettingstarted/>`_

	* For information on using Numpy, Scipy and generally doing scientific computing in Python checkout the `Scipy lecture notes <http://www.scipy-lectures.org/>`_.  The Scipy website also offers as solid introduction to `using Numpy arrays <https://docs.scipy.org/doc/numpy-dev/user/quickstart.html>`_.

	* The `Stackoverflow <http://www.stackoverflow.com>`_ website is an incredible resource for all computing related questions, including simple usage of Python, Scipy and Numpy functions.

	* For users more familiar with Matlab, there is a `Matlab-Numpy cheat sheet <http://mathesaurus.sourceforge.net/matlab-numpy.html>`_ that explains how to translate familiar Matlab commands to Numpy.

===============================================================================
Overview of Data Storage in OpenPNM
===============================================================================

Before creating an OpenPNM simulation it is necessary to give a quick description of how data is stored in OpenPNM; after all, a significant part of OpenPNM is dedicated to data storage and handling.

-------------------------------------------------------------------------------
Python Dictionaries or *dicts*
-------------------------------------------------------------------------------

OpenPNM employs 5 main objects which each store and manage a different type of information or data:

#. **Network**: Manages topological data such as pore spatial locations and pore-to-pore connections
#. **Geometry**: Manages geometrical properties such as pore diameter and throat length
#. **Phase**: Manages thermophysical properties such as temperature and viscosity
#. **Physics**: Manages pore-scale transport parameters such as hydraulic conductance
#. **Algorithm**: Contains algorithms that use the data from other objects to perform simulations, such as diffusion or drainage

We will encounter each of these objects in action before the end of this tutorial.

Each of the above objects is a *subclass* of the Python *dictionary* or *dict*, which is a very general storage container that allows values to be accessed by a name using syntax like:

.. code-block:: python

  	>>> foo = dict()  # Create an empty dict
	>>> foo['bar'] = 1  # Store an integer under the key 'bar'
	>>> foo['bar']  # Retrieve the integer stored in 'bar'
	1

A detailed tutorial on dictionaries `can be found here <http://learnpythonthehardway.org/book/ex39.html>`_.  The *dict* does not offer much functionality aside from basic storage of arbitrary objects, and it is meant to be extended.  OpenPNM extends the *dict* to have functionality specifically suited for dealing with OpenPNM data.  More information about the functionality of OpenPNM's subclassed *dicts* can be found in the :ref:`overall_design`.

-------------------------------------------------------------------------------
*Numpy* Arrays of Pore and Throat Data
-------------------------------------------------------------------------------

All data are stored in arrays which can accessed using standard array syntax.  More details on the data storage scheme are given in :ref:`data_storage`, but the following gives a quick overview:

#. All pore and throat properties are stored in `Numpy arrays <https://docs.scipy.org/doc/numpy-dev/user/quickstart.html>`_.  All data will be automatically converted to a *Numpy* array if necessary.

#. The data for pore *i* (or throat *i*) can be found in element of *i* of an array.  This means that pores and throat have indices which are implied by their position in arrays.  When we speak of retrieving pore locations, it refers to the indices in the *Numpy* arrays.

#. Each property is stored in it's own array, meaning that 'pore diameter' and 'throat volume' are each stored in a separate array.

#. Arrays that store pore data are *Np*-long, while arrays that store throat data are *Nt*-long, where *Np* is the number of pores and *Nt* is the number of throats in the network.

#.  Arrays can be any size in the other dimensions.  For instance, triplets of pore coordinates (i.e. [x, y, z]) can be stored for each pore creating an *Np-by-3* array.

#.  The storage of topological connections is also very nicely accomplished with this 'list-based' format, by creating an array (``'throat.conns'``) that stores which pore indices are found on either end of a throat.  This leads to an *Nt-by-2* array.  The implications and advantages of this storage scheme are discussed further in :ref:`topology`.

-------------------------------------------------------------------------------
OpenPNM Objects: Combining *dicts* and *Numpy* Arrays
-------------------------------------------------------------------------------

OpenPNM objects combine the above two levels of data storage, meaning they are *dicts* that are filled with *Numpy* arrays.  OpenPNM enforces several rules to help maintain data consistency:

#.  When storing arrays in an OpenPNM object, their name (or *dictionary key*) must be prefixed with ``'pore.'`` or ``'throat.'``.

#.  OpenPNM uses the prefix of the *dictionary key* to infer how long the array must be.

#.  The specific property that is stored in each array is indicated by the suffix such as ``'pore.diameter'`` or ``'throat.length'``.

#.  Writing scalar values to OpenPNM objects automatically results in conversion to a full length array filled with the scalar value.

#.  Arrays containing *Boolean* data are treated as *labels*, which are explained later in this tutorial.

The following code snippets give examples of how all these pieces fit together using an **Empty** network as an example:

.. code-block:: python

	>>> import OpenPNM
	>>> import scipy as sp
	>>> net = OpenPNM.Network.Empty(Np=10, Nt=10)  # Instantiate an empty network object with 10 pores and 10 throats
	>>> net['pore.foo'] = sp.ones([net.Np, ])  # Assign an Np-long array of ones
	>>> net['pore.bar'] = range(0, net.Np)  # Assign an Np-long array of increasing ints
	>>> type(net['pore.bar'])  # The Python range iterator is converted to a proper Numpy array
	<class 'numpy.ndarray'>
	>>> net['pore.foo'][4] = 44.0  # Overwrite values in the array
	>>> net['pore.foo'][4]  # Retrieve values from the array
	44.0
	>>> net['pore.foo'][2:6]  # Extract a slice of the array
	array([ 1.,  1., 44.,  1.])
	>>> net['pore.foo'][[2, 4, 6]]  # Extract specific locations
	array([ 1., 44.,  1.])
	>>> net['throat.foo'] = 2  # Assign a scalar
	>>> len(net['throat.foo'])  # The scalar values is converted to an Nt-long array
	10
	>>> net['throat.foo'][4]  # The scalar value was placed into all locations
	2

===============================================================================
Generate a Cubic Network
===============================================================================

Now that we have seen the rough outline of how OpenPNM objects store data, we can begin building a simulation.  Start by importing OpenPNM and the Scipy package:

.. code-block:: python

	>>> import OpenPNM
	>>> import scipy as sp

Next, generate a **Network** by choosing the **Cubic** class, then create an *instance* with the desired parameters:

.. code-block:: python

	>>> pn = OpenPNM.Network.Cubic(shape=[4, 3, 1], spacing=0.0001)

The **Network** object stored in ``pn`` contains pores at the correct spatial positions and connections between the pores according the cubic topology.

* The ``shape`` argument specifies the number of pores in the [X, Y, Z] directions of the cube.  Networks in OpenPNM are always 3D dimensional, meaning that a 2D or "flat" network is still 1 layer of pores "thick" so [X, Y, Z] = [20, 10, 1], thus ``pn`` in this tutorial is 2D which is easier for visualization.

* The ``spacing`` argument controls the center-to-center distance between pores and it can be a scalar or vector (i.e. [0.0001, 0.0002, 0.0003]).

The resulting network looks like:

.. image:: http://i.imgur.com/ScdydO9l.png
   :align: center

This image was creating using `Paraview <http://www.paraview.org>`_, using the instructions given here: `Example in the OpenPNM-Example collection <https://github.com/PMEAL/OpenPNM-Examples/blob/master/IO_and_Visualization/paraview.md>`_

-------------------------------------------------------------------------------
Inspecting Object Properties
-------------------------------------------------------------------------------

OpenPNM objects have additional methods for querying their relevant properties, like the number of pores or throats, which properties have been defined, and so on:

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

More information about these various functions is given in :ref:`overall_design`.  It is also convenient to type ``print(pn)`` at the command line to view a nicely formatted table showing the current state of ``pn``.

-------------------------------------------------------------------------------
Accessing Pores and Throats
-------------------------------------------------------------------------------

One simple but important feature of OpenPNM is the ability to *label* pores and throats.  When a **Cubic** network is created, several labels are automatically created: the pores on each face are labeled 'left', 'right', etc.  These labels can be used as follows:

.. code-block:: python

	>>> pn.pores('left')
	array([0, 3, 6, 9])

The ability to retrieve pore indices is handy for querying pore properties, such as retrieving the pore coordinates of all pores on the 'left' face:

.. code-block:: python

	>>> pn['pore.coords'][pn.pores('left')]
	array([[  5.00000000e-05,   5.00000000e-05,   5.00000000e-05],
	       [  1.50000000e-04,   5.00000000e-05,   5.00000000e-05],
	       [  2.50000000e-04,   5.00000000e-05,   5.00000000e-05],
	       [  3.50000000e-04,   5.00000000e-05,   5.00000000e-05]])

A list of all labels currently assigned to the network can be obtained with:

.. code-block:: python

	>>> pn.labels()
	['pore.all', 'pore.back', 'pore.bottom', 'pore.front', 'pore.internal', 'pore.left', 'pore.right', 'pore.top', 'throat.all']

The existing labels are also listed when an object is printed using ``print(pn)``.  Detailed use of labels is given in :ref:`data_storage`.

===============================================================================
Create a Geometry Object and Assign Geometric Properties to Pores and Throats
===============================================================================

The **Network** ``pn`` does not contain any information about pore and throat sizes at this point.  The next step is to create a **Geometry** object to manage the geometrical properties.

.. code-block:: python

	>>> geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)

This statement contains three arguments:

* ``network`` tells the **Geometry** object which **Network** it is associated with.  There can be multiple networks defined in a given session, so all objects must be associated with a single network.

* ``pores`` and ``throats`` indicate the locations in the **Network** where this **Geometry** object will apply.  In this  tutorial ``geom`` applies to *all* pores and throats, but there are many cases where different regions of the network have different geometrical properties, so OpenPNM allows multiple **Geometry** objects to be created for managing the data in each region, but this is a subject for :ref:`intermediate_usage`.

-------------------------------------------------------------------------------
Add Pore and Throat Size Information
-------------------------------------------------------------------------------

This freshly instantiated **Geometry** object (``geom``) contains no geometric properties as yet.  For this tutorial we'll use the direct assignment of manually calculated values.

We'll start by assigning diameters to each pore from a random distribution, spanning 0 um to 100 um.  The upper limit matches the ``spacing`` of the **Network** which was set to 0.0001 m (i.e. 100 um), so pore diameters exceeding 100 um might overlap with their neighbors.  Using the Scipy ``rand`` function creates an array of random numbers between 0 and 0.0001 that is *Np*-long, meaning each pore is assigned a unique random number

.. code-block:: python

	>>> geom['pore.diameter'] = sp.rand(pn.Np)*0.0001  # Units of meters

We usually want the throat diameters to always be smaller than the two pores which it connects to maintain physical consistency. This requires understanding a little bit about how OpenPNM stores network topology.  Consider the following:

.. code-block:: python

	>>> P12 = pn['throat.conns']  # An Nt x 2 list of pores on the end of each throat
	>>> D12 = geom['pore.diameter'][P12]  # An Nt x 2 list of pore diameters
	>>> Dt = sp.amin(D12, axis=1)  # An Nt x 1 list of the smaller pore from each pair
	>>> geom['throat.diameter'] = Dt

Let's dissect the above lines.

* Firstly, ``P12`` is a direct copy of the **Network's** ``'throat.conns'`` array, which contains the indices of the pore-pair connected by each throat.

* Next, this *Nt-by-2* array is used to index into the ``'pore.diameter'`` array, resulting in another *Nt-by-2* array containing the diameters of the pores on each end of a throat.

* Finally, the Scipy function ``amin`` is used to find the minimum diameter of each pore-pair by specifying the ``axis`` argument as 1, and the resulting *Nt-by-1* array is assigned to ``geom['throat.diameter']``.

* This trick of using ``'throat.conns'`` to index into a pore property array is commonly used in OpenPNM and you should have a second look at the above code to understand it fully.  Refer to :ref:`topology` for a full discussion.

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

The basic geometrical properties of the network are now defined.  The **Geometry** class possesses a method called ``plot_histograms`` that produces a plot of the most pertinent geometrical properties.  The following figure doesn't look very good since the network in this example has only 12 pores, but the utility of the plot for quick inspection is apparent.

.. image:: http://i.imgur.com/xkK1TYfl.png
   :align: center

===============================================================================
Create a Phase Object
===============================================================================

The simulation is now topologically and geometrically defined.  It has pore coordinates, pore and throat sizes and so on.  In order to perform any simulations it is necessary to define a **Phase** object to manage all the thermophysical properties of the fluids in the simulation:

.. code-block:: python

	>>> water = OpenPNM.Phases.GenericPhase(network=pn)

* ``pn`` is passed as an argument because **Phases** must know to which **Network** they belong.

* Note that ``pores`` and ``throats`` are *NOT* specified; this is because **Phases** are mobile and can exist anywhere or everywhere in the domain, so providing specific locations does not make sense.  Algorithms for dynamically determining actual phase distributions are discussed later.

-------------------------------------------------------------------------------
Add Thermophysical Properties
-------------------------------------------------------------------------------

Now it is necessary to fill this **Phase** object with the desired thermophysical properties.  OpenPNM includes a framework for calculating thermophysical properties from models and correlations, but this is covered in :ref:`intermediate_usage`.  For this tutorial, we'll use the basic approach of simply assigning static values as follows:

.. code-block:: python

		>>> water['pore.temperature'] = 298.0
		>>> water['pore.viscosity'] = 0.001

* The above lines utilize the fact that OpenPNM converts scalars to full length arrays, essentially setting the temperature in each pore to 298.0 K.

===============================================================================
Create a Physics Object
===============================================================================

We are still not ready to perform any simulations.  The last step is to define the desired pore-scale physics models, which dictate how the phase and geometrical properties interact to give the *transport parameters*.  A classic example of this is the Hagen-Poiseuille equation for fluid flow through a throat to predict the flow rate as a function of the pressure drop.  The flow rate is proportional to the geometrical size of the throat (radius and length) as well as properties of the fluid (viscosity) and thus combines geometrical and thermophysical properties:

.. code-block:: python

	>>> phys_water = OpenPNM.Physics.GenericPhysics(network=pn, phase=water, geometry=geom)

* As with all objects, the ``Network`` must be specified

* **Physics** objects combine information from a **Phase** (i.e. viscosity) and a **Geometry** (i.e. throat diameter), so each of these must be specified.

* **Physics** objects do not require the specification of which ``pores`` and ``throats`` where they apply, since this information is implied by the ``geometry`` argument which was already assigned to specific locations.

-------------------------------------------------------------------------------
Specify Desired Pore-Scale Transport Parameters
-------------------------------------------------------------------------------

We need to calculate the numerical values representing our chosen pore-scale physics.  To continue with the Hagen-Poiseuille example lets calculate the hydraulic conductance of each throat in the network.  The throat radius and length are easily accessed as:

.. code-block:: python

	>>> R = geom['throat.diameter']/2
	>>> L = geom['throat.length']

The viscosity of the **Phases** was only defined in the pores; however, the hydraulic conductance must be calculated for each throat.  There are several options, but to keep this tutorial simple we'll create a scalar value:

.. code-block:: python

	>>> mu_w = 0.001
	>>> phys_water['throat.hydraulic_conductance'] = 3.14159*R**4/(8*mu_w*L)

Numpy arrays support *vectorization*, so since both ``L`` and ``R`` are arrays of *Nt*-length, their multiplication in this way results in another array that is also *Nt*-long.

===============================================================================
Create an Algorithm Object for Performing a Permeability Simulation
===============================================================================

Finally, it is now possible to run some useful simulations.  The code below estimates the permeability through the network by applying a pressure gradient across and calculating the flux.  This starts by creating a **StokesFlow** algorithm, which is pre-defined in OpenPNM:

.. code-block:: python

	>>> alg = OpenPNM.Algorithms.StokesFlow(network=pn, phase=water)

* Like all the above objects, **Algorithms** must be assigned to a **Network** via the ``network`` argument.

* This algorithm is also associated with a **Phase** object, in this case ``water``, which dictates which pore-scale **Physics** properties to use (recall that ``phys_water`` was associated with ``water``).

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

	>>> Q = alg.rate(pores=pn.pores('front'))
	>>> A = 0.0001*3*1  # Cross-sectional area for flow
	>>> L = 0.0001*4  # Length of flow path
	>>> del_P = 101325  # Specified pressure gradient
	>>> K = Q*mu_w*L/(A*del_P)

The **StokesFlow** class was developed with permeability simulations in mind, so a specific method is available for determining the permeability coefficient that essentially applies the recipe from above.  This method could struggle with non-uniform geometries though, so use with caution:

.. code-block:: python

	>>> K = alg.calc_eff_permeability()

The results (``'pore.pressure'``) are held within the ``alg`` object and must be explicitly returned to the ``air`` object by the user if they wish to use these values in a subsequent calculation.  The point of this data containment is to prevent unintentional overwriting of data.  Each algorithm has a method called ``return_results`` which places the pertinent values back onto the appropriate **Phase** object.

.. code-block:: python

	>>> alg.return_results()

Using Paraview for Visualization, the resulting pressure gradient across the network can be seen:

.. image:: http://i.imgur.com/8aVaH1Sl.png
   :align: center
