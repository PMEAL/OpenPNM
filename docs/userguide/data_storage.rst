.. _network:

###############################################################################
Network Architecture and Data Storage
###############################################################################

OpenPNM utilizes the object oriented capacities of Python by defining a network as an object.  A network object contains both the data that describes the network along with the tools, functions, or methods needed to access this data in ways applicable to the pore network modeling paradigm.  One key feature of this object is that it is completely agnostic about the type of network it describes; a random, cubic or another network topology is stored in exactly the same manner.
The advantages of this are numerous:

1. The OpenPNM framework can be applied to any network or situation
2. All operations performed on the network can be fully generic
3. Only one version of an algorithm or method needs to be developed then can be applied universally
4. The commonality between network types becomes apparent

As the name suggests, pore network modeling borrows significantly from the fields of network and graph theory.  During the development of OpenPNM, it was debated whether existing Python graph theory packages (such as `graph-tool <http://graph-tool.skewed.de/>`_ and `NetworkX <http://networkx.github.io/>`_) should be used to store the network topology.  It was decided that storage of network property data could be handled very efficiently and transparently using simple 1D arrays, which allows for a high degree of code vectorization.  Fortuitously, around the same time as this discussion, Scipy started to include the 'compressed sparse graph' library, which contained numerous graph theory algorithms.  The CSGraph library requires adjacency matrices in a compressed sparse storage scheme, which happens to be how OpenPNM stores network connections as described below.

===============================================================================
Network Architecture
===============================================================================

One of the main design considerations of OpenPNM was to accommodate *all* pore networks (arbitrary dimensionality, connectivity, shape and so on).  Cubic networks are commonly used in pore network modeling, with each pore connected to 6 or 26 neighbors.  This type of network *can* be represented as cubic matrices in numerical simulations, and this has the advantage that it is easily interpreted by human users.  Representing networks this way, however, clearly lacks generality.
Networks extracted from tomographic images, or generated using random pore placements connected by Delaunay tessellations require a different approach.  OpenPNM uses network representation schemes borrowed from graph theory, such as adjacency and incidence matrices, that can be used to represent *all* network topologies.

The basic definitions used by OpenPNM are:

1. Pores can have an arbitrary number of throats, including zero

2. A throat connects exactly two pores, no more and no less

3. Two pores are connected by no more than one throat

A network has a certain number of pores, *Np*, and a certain number of throats, *Nt*.  Typically, *Nt* > *Np* since most pores have more than 1 throat.  If every pore has 1 throat (e.g. the network forms a circular chain), then *Nt* = *Np* - 1.
Of course, *Nt* can be zero but this would not be a useful network.  It can be *unofficially* stated that a network should have at least 2 pores connected by at least 1 throat (*Np* > 1 and *Nt* > 0).

-------------------------------------------------------------------------------
Pore and Throat Numbering
-------------------------------------------------------------------------------

Each pore and throat in the network has a unique ID number.  In OpenPNM the ID number is *implied* by array element location, meaning that any information stored in element *i* of a pore (or throat) property array implicitly applies to pore (or throat) *i*.  Or in other words, finding information about pore (or throat) *i* is accomplished by looking into element *i* of an array.  There is no correspondence between pore number and throat number, meaning that throat *i* may or may not be connected with pore *i*.  Python uses 0-based array indexing so the ID numbers start at 0, which can be a source of confusion when representing connections using sparse representations.  This implicit numbering scheme provides for much easier vectorization of the code, since arrays can be compared with confidence that the elements coincide.
Furthermore, this format also allows for direct indexing based on ID number.  That is, given the ID number(s) for a pore(s), a user can directly lookup a desired property, without having to add a second step of finding where in an array the desired pore information is stored.

.. Note:: The Meaning of the *numbering* property

	One point of confusion that may arise when using OpenPNM is the presence of the pore and throat 'numbering' property.  This property is simply a list of integers between 0 and *Np* or *Nt*.  Its purpose is to facilitate indexing into other property arrays.  For instance, a user can query all pores (or throats) larger than size R and receive a boolean array.  This boolean can then be used to mask the 'numbering' array and produce an array containing a subset of pore (or throat) ID numbers.  The subset can then be passed to algorithms and other methods, such as the inlets to an invasion percolation algorithm.  This approach is a matter of convenience only.  This could be avoided by generating an integer list between 0 and *Np* (or *Nt*) each time the 'numbering' property is used.

-------------------------------------------------------------------------------
'Internal' vs. 'Boundary' Pores
-------------------------------------------------------------------------------

Internal pores and internal throats refer to the throats in which the physical processes occur.  Boundary pores are added to the network to enable numerical calculations that require boundary conditions.  For instance, to simulate diffusion across the network a concentration gradient is created by placing specified concentrations in the boundary pores (Dirichlet conditions).

Boundary pores are not considered part of the physical network; they have no spatial extent thus no volume or length.  They also have no meaningful spatial location, however, for the purposes of visualization they are given coordinates that neighbor the internal pore to which they are connected.  It would be more precise to call them boundary *nodes*, but this leads to other confusions since their properties are stored along with the internal pores in `pore_properties` and `pore_conditions`.

Boundary pores are part of the logical network, thus their ID number and connectivity are vital.  Boundary pores are only connected to internal pores and they are not connected to each other.  Typically, a boundary pore only connects to a single internal pore, but there may be cases where this is not so, such as random networks.  This generally won't impact a simulation.  Internal pores can also be connected to more than one boundary pore.  This can occur when a pore is on an edge or corner of a network and is exposed to multiple boundaries, or can simply result from a confluence of connections, as might occur in a random network.

Throats connecting an internal pore to a boundary pore are considered part of the physical network, so they have spatial extent and location.  The existence of these throats is essential for transmitting the boundary pore information into the physical network.

.. note::

	If surface or boundary effects are of interest, then they must be captured by creating a suitable arrangement of internal pores on the face of the physical domain.  It is not the role of the boundary pores to capture physical processes.

-------------------------------------------------------------------------------
Storing Network Connectivity with Adjacency Matrices
-------------------------------------------------------------------------------
Network topology or connectivity is conveniently and efficiently stored as an `adjacency matrix <http://en.wikipedia.org/wiki/Adjacency_matrix>`_.  An adjacency matrix is a *Np*-by-*Np* 2D matrix.  A non-zero value at location (*i*, *j*) indicates that pores *i* and *j* are connected.  Describing the network in this fashion is one of the main features that allows OpenPNM to be agnostic to the type of network it describes. The fact that all pores and throats are numbered with a unique ID (as described above) is essential to this approach.  Another important feature of the adjacency matrix is that it is highly sparse and can be stored with a variety of sparse storage schemes.  OpenPNM stores the adjacency matrix in the 'COO' or 'IJV' format, which essential stores the coordinates (I,J) and values (V) of the nonzero elements.  Without delving into the details, this approach results in `throat_properties` entry called *'connections'* which is and *Nt*-by-2 array that gives the ID number of the two pores that a given throat connects.  The storage scheme coincides exactly with the storage of all other throat properties.  The details of the OpenPNM implementation of adjacency matrices and other relate issues are given below for the interested reader.

.. Topic:: In Depth: Adjacency and Incidence Matrices

	*Adjacency Matrices*

	When each pore has a unique ID number it is logical to store the network connectivity as a list of the pores to
	which a given pore is connected.  Graph theoreticians have devised an elegant and powerful approach for storing this information, which OpenPNM has adopted, called adjacency matrices.  An adjacency matrix is a sparse 2D matrix of size *Np*-by-*Np*.  A value of 1 is placed at location (*i*, *j*) to indicate that pores *i* and *j* are connected.  In pore networks there is generally no difference between traversing from pore *i* to pore *j* or from pore *j* to pore *i*, so a 1 is also placed at location (*j*, *i*).  This means that determining which pores are connected directly to a given pore (say *i*) can be accomplished by finding the locations of non-zeros in row *i*.  In graph theory terminology this is deemed an *undirected* network, meaning that the *direction* of traversal is immaterial.  The adjacency matrix of an undirected network is symmetric.  Since the adjacency matrix is symmetric it is redundant to store the entire matrix when only the upper (or lower) triangular part is necessary.

	Because pores are generally only connected to nearby pores, the number of throats per pore is a very small fraction of the total number of throats.  This means that there are very few non-zero elements on each row, so the adjacency matrix is highly sparse.  This fact naturally lends itself to sparse storage schemes.  OpenPNM uses uses the IJV sparse storage scheme to store the upper triangular portion of the adjacency matrix.  The *IJV* scheme is simply an *Np*-by-3 array of the (*I*, *J*) coordinates of each non-zero element in the adjacency matrix, along with the corresponding non-zero value (*V*).  (The scipy.sparse module calls this the Coordinate or COO storage scheme, but it is more widely known as IJV).  For example, to denote a value of 1 on row 3 and column 7, the *IJV* storage scheme would include an entry IJV = [3, 7, 1].  Each non-zero element in the adjacency matrix corresponds to a row to the *IJV* array.  Moreover, the number of non-zeros in the upper triangular portion of the adjacency matrix is equal to the number of throats in the network, so the dimensions of the *IJV* array is *Nt*-by-3.  This is not a coincidence; a key feature of the adjacency matrix is that each non-zero element directly corresponds to a throat.  Because throat numbers are implicitly defined by their location in an array, then the IJV sparse storage scheme automatically assigns throat ID numbers when the IJV array is generated.  For instance, when scanning the adjacency matrix from left-to-right, top-to-bottom, the first non-zero element encountered (say at location [0,5]) would be assigned throat number 0, and stored as IJV[0] = [0,5,1].

	One further optimization used by OpenPNM is to drop the V from the IJV format since the non-zeros in the adjacency matrix are all 1.  This results in a *Nt*-by-2 array which is called *connections*.  Any desired throat property array can be appended as a third column to the *connections* array to fully specify the IJV format for use with the scipy.sparse or scipy.csgraph functions.  OpenPNM provides a routine for this operation (``'fill_adjacency_matrix'``), which takes the desired throat property list to insert into *V* as an argument.

	In summary, when storing network connectivity as the upper triangular portion of an adjacency in the IJV sparse storage format, the end result is an *Nt*-by-2 list describing which pores are connected by a given throat.  These connections are a fundamental property associated with each throat in the same way as throat diameter or capillary entry pressure.  This highly distilled storage format minimized memory usage, allows for vectorization of the code, is the most efficient means of generating a sparse matrix, and corresponds perfectly with the storage of other throat properties using the ID number implicitly defined by the list element location.

	*Other Sparse Storage Schemes*

	The IJV storage format corresponds perfectly with the way other throat data is stored in OpenPNM, however some tasks and queries are performed more efficiently using other storage formats.  OpenPNM converts between these formats internally as needed.  For instance, most linear solvers prefer the compressed-sparse-row (CSR) scheme.  Conveniently, the IJV format used by OpenPNM is the fastest way to generate sparse matrices, so conversion, or building of each required sparse format is very efficient.  OpenPNM uses the methods provided by scipy.sparse for these conversions so they are highly optimized and based on C.  OpenPNM contains a method for constructing sparse matrices (called fill_adjacency_matrix) which accepts the storage type as an argument (i.e. 'csr', 'lil', etc).  This method can generate these other formats very quickly since they all derive from the IJV ('coo') format.  For a discussion of sparse storage schemes and the respective merits, see this `Wikipedia article <http://en.wikipedia.org/wiki/Sparse_matrix>`_.

	*Incidence Matrices*

	Another way to represent network connections is an incidence matrix.  This is similar to an adjacency matrix but rather than denoting which pores are connected to which, it denotes which pores are connected to which throats.  An incidence matrix is *Np*-by-*Nt* is size, with *Nt* non-zero elements.  The incidence matrix is useful for quickly querying which throats are connected to a given pore by finding the location of non-zero elements on a row.  Incidence matrices are generated as needed by OpenPNM internally for performing such queries, and the user does not usually interact with them.

===============================================================================
Network Data Storage
===============================================================================
OpenPNM stores two types of information about pores and throats: 'properties' and 'conditions'.  Properties include the geometric and structural aspects of the network, such as pore size and throat length.  Conditions include the thermo-physics and fluids related information such as liquid temperature and gas pressure.  The former information is created by the Geometry modules during network generation, while the latter is produced and altered by the Physics and Algorithm modules.  For instance, an algorithm might calculate the temperature in the network, then a method in the Physics module might use this temperature to calculate temperature dependent liquid viscosity.  There is one important difference between properties and conditions: properties are always vectors of length *Np* for ``pore_properties`` and *Nt* ``throat_properties``, while pore and throat conditions can be either vectors of *Np* and *Nt* respectively, *or* scalars.  The reasons and implications for this will be outlined below.

-------------------------------------------------------------------------------
Pore and Throat Properties
-------------------------------------------------------------------------------
OpenPNM stores all pore and throat properties as Numpy ndarrays.  ndarrays are a numerical data type provided by the Numpy package (which is embedded in the Scipy package) that allow for the type of numerical manipulations that scientists and engineers expect, such as vectorization, slicing, boolean indexing and so on.  Pore properties are stored as arrays of size *Np*-by-*n, where *Np* is the number of pores in the network and *n* is almost always 1, (e.g. pore volume is stored as an *Np*-by-1 array), with a few exceptions (e.g. spatial coordinates are stored as *Np*-by-3 for 3-dimensional space).  Throat properties are almost always stored as *Nt*-by-*m* arrays where *Nt* is the number of throats in the network.  Again, *m* is almost always 1 with a notable exception being the connections property that is discussed in detail above.

As mentioned above, OpenPNM uses implied pore and throat numbering, meaning that the property for pore (or throat) *i* is stored in element *i* of the corresponding property array.

To examine the properties of a network, start by generating a small network of 3-by-3-by-3 as follows:

.. code-block:: python

   import OpenPNM
   pn = OpenPNM.Geometry.Cubic().generate(divisions=[3,3,3],lattice_spacing=[1])

This creates a cubic network with 27 internal pores and 54 internal throats. Additionally, for every 3-by-3 face on the cube, a 3-by-3 set of boundary pores are created with individual boundary throats to corresponding internal pores.  A quick summary of the network data can be displayed as follows:

.. code-block:: python

    print(pn)

	==================================================
	Overview of network properties
	--------------------------------------------------
	Basic properties of the network
	- Number of pores:   81
	- Number of throats: 108

	Pore properties:
		diameter            float64             (81,)
		numbering           int32               (81,)
		volume              float64             (81,)
		seed                float64             (81,)
		coords              float64             (81, 3)
		type                int8                (81,)
	Throat properties:
		volume              float64             (108,)
		diameter            float64             (108,)
		numbering           int32               (108,)
		connections         int32               (108, 2)
		length              float64             (108,)
		seed                float64             (108,)
		type                int8                (108,)

A more detailed description is available with ``pn.print_overview()``.

As can be seen, the default network generation produces several basic pore and throat properties.  Note that the length of the pore and throat property lists correspond to the number of pores and throats in the network (81 and 108 respectively).  Most of the data are stored in 1D arrays, with two exceptions.  The pore property 'coords' gives the spatial location of the pore center in 3D Cartesian coordinates, so each pore requires a set of X, Y and Z values.  The throat property 'connections' gives the ID numbers of the two pores it connects, or in other words it gives the IJ portion of the IJV sparse storage of the adjacency matrix.

These data arrays are stored as part of the network object using Python dictionaries.  A Python dictionary is a form of structured variable where each entry in the dictionary has a { 'key' : <value> } pair.  The 'key' is the name of the of the <value>, and the <value> can be any data type.  In OpenPNM the <values> are all ndarrays.  For example, ``pn.pore_properties['diameter']`` will return the pore diameter values. Similarly, ``pn.throat_properties['diameter']`` returns the throat diameter values.

A quick way to find all properties currently stored in a dictionary is the ``.keys()`` method as follows:

.. code-block:: python

	print(pn.pore_properties.keys())
		['diameter', 'numbering', 'volume', 'seed', 'coords', 'type']

.. note::

	When an ndarray of size *N*-by-1 is used, it is generally preferred to have arrays of shape (N,) rather than (N,1).  There are two reasons for this.  Firstly, in the (N,) form the result of indexing into the array is a scalar, while in the case of (N,1) the result remains a vector and an additional level of index is required to retrieve the actual scalar value.  Secondly, the (N,) case has no transpose so broadcasting during vectorized calculations is failsafe.  In the case of (N,1) there is the possibility of a transposed array of size (1,N) which would lead to an (N,N) result when broadcast.

-------------------------------------------------------------------------------
Pore and Throat Conditions
-------------------------------------------------------------------------------
Pore and throat conditions are very similar to the properties as described above, with one major exception.  'conditions' can be either a vector of length *Np* for pores (and *Nt* for throats), **or** they can be a scalar.  In the case of vector conditions (i.e. one value for each pore or throat) all of the considerations outlined above for 'properties' applies unchanged.  A scalar conditions assumes that this value applies to **all** pores or throats.  For instance, applying a constant temperature to the network can be achieved with:

.. code-block:: python

	pn.pore_conditions['temperature'] = 80.0

Storing this information as a scalar provides significant memory savings by avoiding the redundancy of specifying each pore to have the same temperature.  Fortunately, Numpy is very adapt at 'broadcasting' vectors and scalars together.  This means that a properly vectorized calculation can take a vector or a scalar without any changes to the code.  For instance, to calculate the molar density of the gas in the pores using the ideal gas law, we could write:

.. code-block:: python

	pn.pore_conditions['temperature'] = 80.1
	pn.pore_conditions['pressure'] = 101325
	gas_constant = 8.314
	pn.pore_conditions['molar_density'] = pn.pore_conditions['pressure']/gas_constant/pn.pore_conditions['temperature']

This calculation as shown, with both temperature and pressure as scalars, would produce a scalar value of 'molar_density'.  If, however, either *or* both of 'temperature' and 'pressure' were vectors (i.e. a value for each pore), then the 'molar_density' would be calculated in *exactly* the same way, but the result would be a vector.  The only caveat is that all vectors involved must be the same length.

.. Topic:: Upcoming Feature

	**Special Features of the OpenPNM Dictionaries**

	The dictionaries used in OpenPNM have been sub-classed from the general Python implementation.  Since so many operations in OpenPNM depend on vectorized code, it is imperative that all ``pore_properties`` arrays are a consistent length (and similarly for ``throat_properties``).  Pythons native dictionary class has been extended to include a check for array shape prior to adding or overwriting arrays.  The *self-protecting* properties of this dictionary will be expanded in future releases as the develops.

The ``pore_conditions`` and ``throat_conditions`` arrays are also written in dictionaries, but as mentioned above, scalar values are allowed.  The dictionary class in OpenPNM allows this, as well as allowing a scalar to be expanded to an *Np* or *Nt* vector.  It will not allow vectors of lengths other than these.

-------------------------------------------------------------------------------
Mandatory Pore and Throat Properties
-------------------------------------------------------------------------------
The default behavior of the GenericGeometry generator produces several pore and throat properties based on commonly used assumptions.  Only a few of these properties are truly essential to defining the pore network.

**'connections' and 'coords'**

The spatial position of each pore is obviously a defining feature of a given pore network, so the 'coords' pore property is essential.  Equally essential to defining a network is the 'connections' throat property since this describes how the pores are connected or networked.  From a physical point of view, these are the only properties required to define a basic (though not very functional) network.  With this information it would be possible to generate a 3D images of the pore and throat network.

**'type' and 'numbering'**

The 'type' and 'numbering' properties are also considered mandatory since OpenPNM relies on these for various internal calculations and network queries.

The 'numbering' array is actually somewhat redundant since pore and throat numbers are implicitly defined by their array location.  This array is quite useful for boolean mask logic to find pores that meet a specific criteria.  For instance, to find all pores whose positions are below the mean:

.. code-block:: python

	dia_mean = sp.mean(pn.pore_properties['diameter'])
	mask = pn.pore_properties['diameter'] < dia_mean
	small_pores = pn.pore_properties['numbering'][mask]
	print(small_pores)

	[ 0  1  2  3  4  5  6  7  8 27 28 29 36 37 38 45 46 47 54 55 56 63 64 65 66 67 68 69 70 71]

The 'type' property is used by OpenPNM to differentiate between internal and boundary pores (and throats).  A 'type' value of zero indicates an internal pore, and a value > 0 indicates a boundary pore.  Boundary pores are further distinguished by values between 1 and 6 to indicate on which boundary they lie: 1 and 6 for z-faces, 2 & 5 for x-faces and 3 & 4 for y-faces.  This convention was inspired by the number on dice, where opposite sides all add up to 7.  Obviously, this numbering boundary pores in this way implies a cubic network domain, which may not always be the case. For example, let's determine on which faces the low_pores reside:

.. code-block:: python

	>>> print(pn.pore_properties['type'][low_pores])
	[0 0 0 0 0 0 0 0 0 2 2 2 5 5 5 3 3 3 4 4 4 1 1 1 1 1 1 1 1 1]

Throats are by definition always internal to the network, but they also have a 'type' property.  If throats are connected to a boundary pore, then they adopt this pores type, otherwise they are 0.
-------------------------------------------------------------------------------
Common Pore and Throat Properties
-------------------------------------------------------------------------------
The GenericGeometry class includes several methods that produce some additional pore and throat properties beyond the mandatory ones described above.  These including this like 'diameter' and 'volume'.  The docstrings for the methods in the GenericGenerator are provided below, with small blurbs about what properties are created at each step and how.

-------------------------------------------------------------------------------
Adding New Pore and Throat Dictionary Entries
-------------------------------------------------------------------------------
Adding a new entry into either of the *properties* or *conditions* dictionaries is very straight-forward.  For instance, creating a throat property called 'aspect_ratio' is as simple as:

.. code-block:: python

	Nt = pn.get_num_throats()
	values = sp.random.rand(Nt,)*4 + 1 # 1 < ratios < 5
	pn.throat_properties['aspect_ratio'] = values

The length of the array generated here is *Nt*, so an aspect ratio is assigned to each throat.
Attempts to add entries of the wrong size would be intercepted by the dictionary class to prevent corruption of the network data.

===============================================================================
Querying Network Data and Properties
===============================================================================

The OpenPNM network object not only stores the network data, but also contains numerous methods for extracting information about the network from that data.
A full listing of the available methods can be found in the :ref:`Network function reference <network_ref>`.

Two particularly useful methods are get_neighbor_pores and get_neighbor_throats.
Their functionality and applicability will be outlined here, which will also give insight to the usage of the other related methods.

In it's basic form `get_neighbor_pores(i)` returns the ID number of the pores directly connected to pore i.
This is useful for checking the status of a given pores neighbors


This function however also takes a list of pores and returns a cumulative list of the neighbors


.. code-block:: python

	>>> print('There are',pn.get_num_pores(),'pores in the network.')
	There are 81 pores in the network.
	>>> print('There are',pn.get_num_throats(),'throats in the network.')
	There are 108 throats in the network.
	>>> print('Pore 5 has the following neighbors:',pn.get_neighbor_pores(5))
	Pore 5 has the following neighbors: [ 2  4  8 14 37 68]
	>>> print('Pore 5 has the following throats:',pn.get_neighbor_throats(5))
	Pore 5 has the following throats: [ 6 11 14 15 64 95]
	>>> print('Pore 5 has',pn.get_num_neighbors(5),'neighbors.')
	Pore 5 has [6] neighbors.
	>>> print('Throat 6 connects pore',pn.get_connected_pores(6)[0],'to pore',pn.get_connected_pores(6)[1],'.')
	Throat 6 connects pore 2 to pore 5
	>>> print('Throat',pn.get_connecting_throat(0,1),'connects pore 0 to pore 1.')
	Throat [0] connects pore 0 to pore 1.


