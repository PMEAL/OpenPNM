. _geometry:

###############################################################################
Pore and Throat Geometry
###############################################################################
In OpenPNM the pore and throat geometry are defined separately from the Network topology.  In other words, creating a network simply places pores at certain coordinates and connects them in a certain pattern.  It is the job of the Geometry object to calculate the physical properties of the pores and throats (i.e. sizes, volumes, lengths, etc), based on a specific pore or throat model (i.e. sphere, cuboid, cylinder, etc).  

Geometry objects are `custom built` by the user in 2 ways: Firstly, the user can select which physical models they wish to use and how the physical properties of should be calculated.  Secondly, the user can add their own physical models and property calculation methods to the list supplied with OpenPNM to extend the framework.

===============================================================================
Generating Geometry
===============================================================================
The most general way to generate a Geometry object is as follows:

.. code::

    >> pn = OpenPNM.Network.TestNet()  # This generates a 5x5x5 cubic network for testing purposes
    >> geom = OpenPNM.Geometry.GenericNetwork(network=pn, name='custom_1', locations='all')
	
There are 3 arguments sent to GenericGeometry here.  Firstly, the Geometry object must be associated with a network (``network=pn``).  This gives the Geometry object access to the network topology such as the number of pores and throats, where they are located in space, and how they are connected.  This is required for something like throat length which depends on the distance between two neighboring pores.  Secondly, each Geometry object (and all objects in OpenPNM) must be given a unique name (``name='custom_1'``).  This makes is possible for a human to differentiate between multiple different Geometry objects by simply checking their ``name`` attribute.  Finally, the locations argument (``locations='all'``) tells the Geometry object to which pores it applies.  A single Network can have a many different Geometry objects associated with it.  For instance a region of low permeability might be embedded in the middle of the domain, so the Geometry object for this region would calculate much smaller pores.  In this case the locations were set by assigning ``geom`` to any pore with the label ``'all'``` which by definition is all pores.  Pore and throat labels are discussed in the :ref:`introduction`.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Building a Custom Geometry Object
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Once an empty GenericGeometry object has been instantiated, the next step is to add methods to the object that calculate the appropriate pore and throat properties.  GenericGeometry contains the ``add_method`` function for this purpose.  It is typical to assign a random seed to each pore in the network which is subsequently used to calculate pore diameters from a statistical distribution.  OpenPNM.Geometry comes with a submodule called ``pore_seed`` that contains several methods that can be used.  The desired method is added to the Geometry object as follows:

.. code::

    >> geom.add_method(prop='pore_seed')











