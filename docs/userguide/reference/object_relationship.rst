.. _object_relationship:

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Class Hierarchy
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The design of OpenPNM is to separate different types of properties between different objects.  There are 5 types: **Network**, **Geometry**, **Phase**, **Physics**, and **Algorithms**.  Each of these are described in more detail below, but their names clearly indicate what sort of data or calculations are assigned to each.

===============================================================================
Core
===============================================================================
Using the Object Oriented Programming (OOP) approach, each of the above objects actually descends from a common class called **Core**.  The **Core** class defines the majority of the functionality, which is then enhanced and extended by each descendent.

**Core** is a subclass of the Python Dictionary or ``dict``.  A ``dict`` is a very handy data structure that can store any piece of data by name, using the following:

.. code-block:: python

    >>> # Instantiate a dict and add some values by name
    >>> foo = dict()
    >>> foo['an_int'] = 1
    >>> foo['a_list'] = [1, 2, 3]
    >>> foo['a_string'] = 'bar'
    >>> # And data can be retrieved by name
    >>> foo['an_int']
    1
    >>> foo['a_list']
    [1, 2, 3]

They Python ``dict`` class comes with a variety of methods for adding, removing, and inspecting the data stored within.  The following command will generate a list of all these methods, which include things like ``pop`` for removing items from the dictionary, and ``keys`` for listing all the current dictionary entries.

.. code-block:: python

    >>> methods = [item for item in dir(foo}) if not item.startswith('_')]

The **Core** class possess all of these methods, plus another few dozen methods that were added by OpenPNM.  These additional methods also pertain to the manipulation of data, but are specific to the types of data used in OpenPNM.

1.  ``props`` and ``labels``
Returns a list of which properties or labels exist in the dictionary.  These methods are basically the same as the ``keys`` method, but return a subset of the entries.  Any arrays of Boolean type are considered labels, while all other are properties.  The returned list outputs a nicely formatted table to the command line when it is printed.

2.  ``num_pores`` and ``num_throats`` (plus ``Np`` and ``Nt``)
Returns the number of pores or throats that the object controls.  Both optionally accept a list of labels and returns the number of pores or throats possessing those labels.  There is a ``mode`` argument which allows control over how the label query is performed.  ``Np`` and ``Nt`` are short-cuts that return the total number of pores or throats.

3.  ``pores`` and ``throats`` (plus ``Ps`` and ``Ts``)
Returns a list of pore or throat indices.  Both optionally accept a list of labels and returns only a list of pores or throats possessing those labels.  There is a ``mode`` argument which allows control over how the label query is performed.  ``Ps`` and ``Ts`` are short-cuts that return ALL of the pore or throat indices.

4.  ``tomask`` and ``toindices``
These methods allow the conversion between numeric indices and Boolean masks.

5.  ``map_pores`` and ``map_throats`` (plus ``Pnet`` and ``Tnet``)
Each **Core** object has it's own internal numbering scheme, so these methods are for converting the pore or throat indices from one object to another.  Practically speaking this usually means mapping from a **Geometry** or **Physics** object onto the **Network**, so ``Pnet`` and ``Tnet`` are short-cuts for retrieving a list of pore or throat indices on the network.

6.  ``network``, ``geometries``, ``phases``, ``physics``, and ``controller``
When each object is instantiated it is associated with the other objects within the simulation.  These methods allow for retrieval of these other objects.

7.  ``interpolate_data``
Data is often calculated or assigned to pores or throats only.  This method enables the conversion of data between these.

8.  ``check_data_health``
Checks whether any data on the object is not well formed, such as containing NaNs, or infs.  This is handy be running an algorithm to ensure that all necessary properties have been defined everywhere.

9.  ``add_model``, ``regenerate``, and ``models``
The ``models`` attribute actually contains a nested dictionary which stores all the information related to the pore-scale models.  This is described elsewhere in detail.  ``add_model`` and ``regenerate`` are wrapper or helper methods to provide quicker access to the ``add`` and ``regenerate`` methods of the ``models`` dict.

10.  ``name``
Contains a unique string identifier for the object.  It can be specified or assigned at will, but no to objects can have the same name.

===============================================================================
Network
===============================================================================
1.  ``check_geometry_health``
Inspects that all pores and throats have been assigned to a **Geometry** object.

2.  ``check_network_health``
Performs a suite of topological checks for ill conditioned networks (disconnected pores, duplicate throats, etc.)

3.  ``clone_pores``, ``connect_pores``, ``extend``, ``stitch``, ``trim``
These are topological manipulation methods that are used to add or remove pores and throats from the network.  These are helper methods for the actual functions in **Network.tools**.

4.  ``find_neighbor_pores``, ``find_neighbor_throats``, ``find_nearby_pores``, ``find_connected_pores`` and ``find_connecting_throat``
These methods can be used to query the neighborhood around a given set of pores.

5.  ``create_adjacency_matrix``, ``create_incidence_matrix``
Returns a *Scipy Sparse* array describing the topology of the network.

6.  ``find_clusters`` and ``find_clusters2``
Finds connected clusters of pores based on a given list of Boolean values.  The 2nd generation of this algorithm has more options that the original, which was kept for backwards compatibility.

7.  ``domain_area``, ``domain_length``, ``domain_bulk_volume``, ``domain_pore_volume``
These calculate the bulk dimensions of the domain.

===============================================================================
Geometry
===============================================================================
1.  ``set_locations``
When instantiating a **Geometry** object it is normal to specify which pores and throats it applies to.  These can be adjusted after the fact with this method.

===============================================================================
Phase
===============================================================================
1.  ``check_physics_health``
Inspects that all pores and throats have been assigned to a **Physics** object.

2.  ``check_mixture_health``
Mixtures are not fully implemented yet, but this makes sure all mole fractions sum to 1.

===============================================================================
Physics
===============================================================================
1.  ``set_locations``
When instantiating a **Physics** object it is normal to specify which pores and throats it applies to.  These can be adjusted after the fact with this method.

2.  ``parent_phase``
``phases`` gives the ability to find a list of all **Phases** in the simulation, but this method returns a handle to the specific **Phase** it's associated with.

===============================================================================
Algorithms
===============================================================================
Depending on the **Algorithm** in question, the additional methods can vary.  Most have:

1.  ``setup``
This method is called to specify some of the optional parameters

2.  ``set_boundary_conditions``
Used to specify the boundary conditions of the simulation.  Some methods also include ``set_inlets`` and ``set_outlets``.  
