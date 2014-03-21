.. _geometry:

###############################################################################
Pore and Throat Geometry
###############################################################################
In OpenPNM the pore and throat geometry are defined separately from the **Network** topology.  In other words, creating a network simply places pores at certain coordinates and connects them in a certain pattern.  It is the job of the **Geometry** object(s) to calculate the physical properties of the pores and throats (i.e. sizes, volumes, lengths, etc), based on a specific pore or throat model (i.e. sphere, cuboid, cylinder, etc).  

.. note:: 

	Fluid, Geometry and Physics modules are designed to function identically, so once you're familiar with the usage of one then all the others should be similar.  

===============================================================================
What is a Geometry Object?
===============================================================================

**Geometry** objects have one main function in OpenPNM.  They contain the models the user wishes to use to calculate pore and throat properties.  

1. The user selects which physical models they wish to use and how the physical properties of should be calculated.  

2. Secondly, the user can add their own physical models and property calculation methods to the list supplied with OpenPNM to extend the framework. 

===============================================================================
Generating Geometry
===============================================================================
The most general way to generate a **Geometry** object is as follows:

.. code-block:: python

	pn = OpenPNM.Network.TestNet()  # This generates a 5x5x5 cubic network for testing purposes
	geom = OpenPNM.Geometry.GenericNetwork(network=pn, name='custom_1')
	
There are 2 arguments sent to ``GenericGeometry`` here.  Firstly, the **Geometry** object is associated with a **Network** with ``network=pn``.  This gives the **Geometry** object access to the network topology such as the number of pores and throats, where they are located in space, and how they are connected.  This is required for something like throat length which depends on the distance between two neighboring pores.  Secondly, each **Geometry** object (and all objects in OpenPNM) must be given a unique name (``name='custom_1'``).  This makes is possible for a human to differentiate between multiple different **Geometry** objects by simply checking their ``name`` attribute.  

In the above example, the **Geometry** object ``geom`` will apply to all pores and throats.  A single **Network** can have a many different **Geometry** objects associated with it.  For instance a region of low permeability might be embedded in the middle of the domain, so the **Geometry** object for this region would calculate much smaller pores.  In this case it is necessary to initialize each **Geometry** object with a list of which pores and throats it applies to.  Assuming that pores and throats for two domains have already been given labels of 'subdomain1' and 'subdomain2', the following procedure would generate two **GeometryGG objects and apply them to the correct locations.  

.. code-block:: python

	pn = OpenPNM.Network.TestNet()  # This generates a 5x5x5 cubic network for testing purposes
	geom = OpenPNM.Geometry.GenericNetwork(network=pn, name='custom_1')

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Building a Geometry Object
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Once a **Geometry** object has been instantiated, the next step is to add methods to the object that calculate the appropriate pore and throat properties.  The ``GenericGeometry`` class contains the ``add_method` function for this purpose.  It is typical to assign a random seed to each pore in the network which is subsequently used to calculate pore diameters from a statistical distribution.  The **Geometry** module comes with a submodule called **pore_seed** that contains several methods that can be used.  The desired method is added to the **Geometry** object as follows:

.. code-block:: python

	geom.add_method(prop='pore_seed',method='random')
	geom.add_method(prop='pore_seed',propname='second_seed',method='random')
	
The first line above looks into the **pore_seed** submodule and finds a method named ``random``.  It attaches this method to itself under the attribute name ``pore_seed`` because the default is to use the name of ``prop``.  The **Geometry** object now knows how to generate pore seed values when they are needed.  The second line above does precisely the same procedure, except it stores the method under the attribute name ``second_seed``.  The ``propname`` argument makes it possible to store as many methods on the object as desired without overwriting previously attached methods.  So ``prop`` is the name of the submodule where ``add_method`` looks to find the method, but ``propname`` is where the method is stored.  

A critical feature of this method addition approach is that none of the arguments to the method can be changed after it has been added to the object.  This is done deliberately using Python's 'partial function' feature.  

In the above case the seed generated by each method will differ since the state of the random number generator was not specifically.  Many methods, including the ``random`` mothod in **pore_seed** accept or require additional parameters. In the case of ``random`` it is possible to send a seed value which initializes Scipy's random number generator to the specified state as follows:

.. code::

	geom.add_method(prop='pore_seed',method='random',seed=10)

Attaching the ``pore_seed`` method to the **Geometry** object in this way will always result in the same random numbers being placed inside each pore since the generator will always be initiated with to the same starting state.  

The above procedure is repeated for all the desired methods.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Generating or Regenerating Geometry Data
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Once the **Geometry** object has been built and contains all the desired property models it is necessary to actually run all these methods to calculate their results.   When the time comes to generate the pore and throat size data (or regenerate it) one *can* manually run each method as:

.. code-block:: python

	geom.pore_seed()

If a large number of method have been added and/or they have custom names this can be an annoying task.  To avoid this, each time ``add_method`` is called it appends the 'propname' to a private list of attached methods.  The ``GenericGeometry`` class includes a method called ``regenerate`` which simply scans through this list and calls each method.  The items in the list are stored in the order they were called in, and the methods are invoked in that order.  It is possible to regenerate only some methods by sending their attribute name ('propname') to the ``regenerate`` method as a list of strings.  It is also possible to exclude certain method from being run listing them in the ``exclude`` argument, if for some reason you don't want to regenerate certain properties.  

===============================================================================
Available Property Estimation Models
===============================================================================

For a complete list of available pore scale geometry models see the :ref:`Function Reference <geometry_ref>`.
