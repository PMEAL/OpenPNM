.. _overall_design:

###############################################################################
Overall Design
###############################################################################

The design of OpenPNM is to separate different types of properties between different objects.  There are 5 types: **Network**, **Geometry**, **Phase**, **Physics**, and **Algorithms**.  Each of these are described in more detail below, but their names clearly indicate what sort of data or calculations are assigned to each.

The image below outlines how each of the main objects in OpenPNM descend from the Python ``dict`` class.  Using the Object Oriented Programming (OOP) paradigm, each of the main OpenPNM objects actually descend from a common class called **Core**.  The **Core** class defines the majority of the functionality, which is then enhanced and extended by each descendent.  This extra functionality is explored in more detail for each main object in the sections below.

.. image:: http://i.imgur.com/k4Far1W.png
   :width: 800 px
   :align: center

===============================================================================
Core
===============================================================================

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

The Python ``dict`` class comes with a variety of methods for adding, removing, and inspecting the data stored within.  The following command will generate a list of all these methods, which include things like ``pop`` for removing items from the dictionary, and ``keys`` for listing all the current dictionary entries.   The following line inspects the ``dict`` object and retrieves a list of all methods it has available to it, which is stored in ``methods``.

.. code-block:: python

    >>> methods = [item for item in dir(foo) if not item.startswith('_')]

The **Core** class possess all of these methods, plus another few dozen methods that were added by OpenPNM.  These additional methods also pertain to the manipulation of data, but are specific to the types of data used in OpenPNM.

-------------------------------------------------------------------------------
1.  Querying Defined Properties and Labels
-------------------------------------------------------------------------------
**Returns a list of which properties or labels exist in the dictionary**

These methods are basically the same as the ``keys`` method, but return a subset of the entries.  Any arrays of Boolean type are considered labels, while all other are properties.  The returned list outputs a nicely formatted table to the command line when it is printed.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Base.Core.props()
    ~OpenPNM.Base.Core.labels()

-------------------------------------------------------------------------------
2. Counting Pores and Throats
-------------------------------------------------------------------------------
**Returns the number of pores or throats that the object controls**

Both optionally accept a list of labels and returns the number of pores or throats possessing those labels.  There is a ``mode`` argument which allows control over how the label query is performed.  ``Np`` and ``Nt`` are short-cuts that return the total number of pores or throats.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Base.Core.num_pores()
    ~OpenPNM.Base.Core.num_throats()
    ~OpenPNM.Base.Core.Np
    ~OpenPNM.Base.Core.Nt

-------------------------------------------------------------------------------
3.  Retrieving a List of Specific Pores and Throats
-------------------------------------------------------------------------------
Returns a list of pore or throat indices.  Both optionally accept a list of labels and returns only a list of pores or throats possessing those labels.  There is a ``mode`` argument which allows control over how the label query is performed.  ``Ps`` and ``Ts`` are short-cuts that return ALL of the pore or throat indices.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Base.Core.pores()
    ~OpenPNM.Base.Core.throats()
    ~OpenPNM.Base.Core.Ps
    ~OpenPNM.Base.Core.Ts

-------------------------------------------------------------------------------
4.  Converting Between Masks and Indices
-------------------------------------------------------------------------------
These methods allow the conversion between numeric indices and Boolean masks.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Base.Core.tomask()
    ~OpenPNM.Base.Core.toindices()

-------------------------------------------------------------------------------
5.  Mapping Pore and Throat Indices Between Objects
-------------------------------------------------------------------------------
Each **Core** object has it's own internal numbering scheme, so these methods are for converting the pore or throat indices from one object to another.  Practically speaking this usually means mapping from a **Geometry** or **Physics** object onto the **Network**, so ``Pnet`` and ``Tnet`` are short-cuts for retrieving a list of pore or throat indices on the network.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Base.Core.map_pores()
    ~OpenPNM.Base.Core.map_throats()
    ~OpenPNM.Base.Core.Pnet
    ~OpenPNM.Base.Core.Tnet

-------------------------------------------------------------------------------
6.  Looking Up Other Objects in the Simulation
-------------------------------------------------------------------------------
When each object is instantiated it is associated with the other objects within the simulation.  These methods allow for retrieval of these other objects.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Base.Core.network
    ~OpenPNM.Base.Core.geometries
    ~OpenPNM.Base.Core.phases
    ~OpenPNM.Base.Core.physics

-------------------------------------------------------------------------------
7.  Interpolating Between Pore and Throat Data
-------------------------------------------------------------------------------
Data is often calculated or assigned to pores or throats only.  This method enables the conversion of data between these.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Base.Core.interpolate_data()

-------------------------------------------------------------------------------
8.  Check the Health of all Data Arrays
-------------------------------------------------------------------------------
Checks whether any data on the object is not well formed, such as containing NaNs, or infs.  This is handy be running an algorithm to ensure that all necessary properties have been defined everywhere.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Base.Core.check_data_health()

-------------------------------------------------------------------------------
9.  Using Pore-Scale Models
-------------------------------------------------------------------------------
The ``models`` attribute actually contains a nested dictionary which stores all the information related to the pore-scale models.  This is described elsewhere in detail.  ``add_model`` and ``regenerate`` are wrapper or helper methods to provide quicker access to the ``add`` and ``regenerate`` methods of the ``models`` dict.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Base.Core.add_model()
    ~OpenPNM.Base.Core.regenerate()

-------------------------------------------------------------------------------
10.  Find and Set the Object's Name
-------------------------------------------------------------------------------
Contains a unique string identifier for the object.  It can be specified or assigned at will, but no to objects can have the same name.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Base.Core.name

===============================================================================
Network
===============================================================================

-------------------------------------------------------------------------------
1.  Check the Health of Associated Geometry Objects
-------------------------------------------------------------------------------
Inspects that all pores and throats have been assigned to a **Geometry** object.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Network.GenericNetwork.check_geometry_health()

-------------------------------------------------------------------------------
2.  Check the Health of the Network Topology
-------------------------------------------------------------------------------
Performs a suite of topological checks for ill conditioned networks (disconnected pores, duplicate throats, etc.)

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Network.GenericNetwork.check_network_health()

-------------------------------------------------------------------------------
3.  Manipulate Pore Topology
-------------------------------------------------------------------------------
These are topological manipulation methods that are used to add or remove pores and throats from the network.  These are helper methods for the actual functions in **Network.tools**.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Network.GenericNetwork.clone_pores()
    ~OpenPNM.Network.GenericNetwork.connect_pores()
    ~OpenPNM.Network.GenericNetwork.extend()
    ~OpenPNM.Network.GenericNetwork.stitch()
    ~OpenPNM.Network.GenericNetwork.trim()

-------------------------------------------------------------------------------
4.  Query Neighborhood
-------------------------------------------------------------------------------
These methods can be used to query the neighborhood around a given set of pores.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Network.GenericNetwork.find_neighbor_pores()
    ~OpenPNM.Network.GenericNetwork.find_neighbor_throats()
    ~OpenPNM.Network.GenericNetwork.find_nearby_pores()
    ~OpenPNM.Network.GenericNetwork.find_connected_pores()
    ~OpenPNM.Network.GenericNetwork.find_connecting_throat()

-------------------------------------------------------------------------------
5.  Adjacency and Incidence Matrices
-------------------------------------------------------------------------------
Returns a *Scipy Sparse* array describing the topology of the network.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Network.GenericNetwork.create_adjacency_matrix()
    ~OpenPNM.Network.GenericNetwork.create_incidence_matrix()

-------------------------------------------------------------------------------
6.  Search for Clusters of Pores
-------------------------------------------------------------------------------
Finds connected clusters of pores based on a given list of Boolean values.  The 2nd generation of this algorithm has more options that the original, which was kept for backwards compatibility.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Network.GenericNetwork.find_clusters()
    ~OpenPNM.Network.GenericNetwork.find_clusters2()

-------------------------------------------------------------------------------
7.  Query the Domain Size
-------------------------------------------------------------------------------
These calculate the bulk dimensions of the domain.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Network.GenericNetwork.domain_area()
    ~OpenPNM.Network.GenericNetwork.domain_length()
    ~OpenPNM.Network.GenericNetwork.domain_bulk_volume()
    ~OpenPNM.Network.GenericNetwork.domain_pore_volume()

===============================================================================
Geometry
===============================================================================

-------------------------------------------------------------------------------
1.  Assign Geometry to Specific Pores and Throats
-------------------------------------------------------------------------------
When instantiating a **Geometry** object it is normal to specify which pores and throats it applies to.  These can be adjusted after the fact with this method.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Geometry.GenericGeometry.set_locations()

===============================================================================
Phase
===============================================================================

-------------------------------------------------------------------------------
1.  Check the Health of Associated Physics Objects
-------------------------------------------------------------------------------
Inspects that all pores and throats have been assigned to a **Physics** object.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Phases.GenericPhase.check_physics_health()

-------------------------------------------------------------------------------
2.  Check the Health of a Mixture Phase
-------------------------------------------------------------------------------
Mixtures are not fully implemented yet, but this makes sure all mole fractions sum to 1.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Phases.GenericPhase.check_mixture_health()

===============================================================================
Physics
===============================================================================

-------------------------------------------------------------------------------
1.  Assign Physics to Specific Pores and Throats
-------------------------------------------------------------------------------
When instantiating a **Physics** object it is normal to specify which pores and throats it applies to.  These can be adjusted after the fact with this method.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Physics.GenericPhysics.set_locations()

-------------------------------------------------------------------------------
2.  Lookup the Parent Phase
-------------------------------------------------------------------------------
The ``phases`` method of the **Core** class gives the ability to find a list of all **Phases** in the simulation, but this method returns a handle to the specific **Phase** it's associated with.

.. autosummary::
    :toctree: generated/

    ~OpenPNM.Physics.GenericPhysics.parent_phase()

===============================================================================
Algorithms
===============================================================================

Depending on the **Algorithm** in question, the additional methods can vary.  Most have:

-------------------------------------------------------------------------------
1.  Specifying Setup Parameters
-------------------------------------------------------------------------------
This method is called to specify some of the optional parameters

-------------------------------------------------------------------------------
2.  Setting Boundary Conditions
-------------------------------------------------------------------------------
Used to specify the boundary conditions of the simulation.  Some methods also include ``set_inlets`` and ``set_outlets``.
