.. _network:

===============================================================================
Network
===============================================================================
This module is the heart of OpenPNM.  It contains the ``GenericNetwork`` class which possesses a suite of network query methods, based the graph theory concepts of adjacency and incidence matrices.  The methods in ``GenericNetwork`` are fully agnostic to the type and topology of network due the generalized way that OpenPNM stores data.  This is explained in more detail in the :ref:`here<data_storage>`. 

All OpenPNM objects descend from the Python ``dict`` object, and have a number of additional methods for dealing with the OpenPNM related data.  The GenericNetwork class is special in that it has a large number of *additional* methods to work with the topological data. 

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Basic Usage
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The ``GenericNetwork`` class on its own has no topology.  If you instantiate a ``GenericNetwork`` it will have no pores or throats.  You can get a quick overview of the network properties by 'printing' it on the command line:

>>> print(pn)
------------------------------------------------------------
OpenPNM.Network.GenericNetwork: 	GenericNetwork_GnSpz
------------------------------------------------------------
#     Properties                          Valid Values
------------------------------------------------------------
1     pore.coords                             0 / 0    
2     throat.conns                            0 / 0    
------------------------------------------------------------
#     Labels                              Assigned Locations
------------------------------------------------------------
1     pore.all                            0         
2     throat.all                          0         
------------------------------------------------------------

As can be seen, a basic empty network has 0 pore coordinates and 0 throat connections, and the label 'all' exists but is applied nowhere.  It is also possible to just get the number of pores or throats on the object:

>>> pn = OpenPNM.Network.GenericNetwork()
>>> pn.num_pores()
0
>>> pn.num_throats()
0

The Network module contains numerous subclasses of ``GenericNetwork``, which possess the code for actually generating specific network topologies (e.g. cubic, random, etc).  All subclasses derive from ``GenericNetwork`` so have its methods, as well as any additional methods relevant to the specific topology.  Generating a standard Cubic network is accomplished with:

>>> pn = OpenPNM.Network.Cubic(shape=[3,3,3],name='demo')
>>> print(pn)
------------------------------------------------------------
OpenPNM.Network.Cubic: 	demo
------------------------------------------------------------
#     Properties                          Valid Values
------------------------------------------------------------
1     pore.coords                            27 / 27   
2     pore.index                             27 / 27   
3     throat.conns                           54 / 54   
------------------------------------------------------------
#     Labels                              Assigned Locations
------------------------------------------------------------
1     pore.all                            27        
2     pore.back                           9         
3     pore.bottom                         9         
4     pore.front                          9         
5     pore.internal                       27        
6     pore.left                           9         
7     pore.right                          9         
8     pore.top                            9         
9     throat.all                          54        
------------------------------------------------------------

The print-out of the network information shows that it has 27 pores and 54 throats, with properties of 'pore.coords', 'pore.index' and 'throat.conns'.  Because the ``Cubic`` class only generates the topology, there is not any information about pore and throat sizes.  The only requirements of a topology are that the pores have spatial locations (given by 'pore.coords') and throats know which two pores they connect ('throat.conns').  ('pore.index' is required for other purposes).  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Topology Queries
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The OpenPNM subclass of the Python ``dict`` has numerous additional methods that are all available to all the main OpenPNM objects.  The GenericNetwork class has an additional suite of methods that are specifically relating to querying the Network topology, such as finding the neighbors of a pore, or finding the throat that connects 2 pores:

>>> pn.find_neighbor_pores(pores=[0])
array([1, 3, 9])
>>> pn.find_connecting_throat(P1=[0,0,0],P2=[1,3,9])
[[0], [18], [36]]
>>> pn.find_connected_pores(throats=[0,18,36])
array([[0, 1],
       [0, 3],
       [0, 9]])

The best way to explore the available methods is to use an IDE or editor that support the autocomplete function, such as Spyder.  This way, you can type ``pn.`` and a pop-up list of available methods will appear.  Extensive documentation is also included inside the OpenPNM code itself in the form of 'docstrings' which will be interpreted by Spyder and shown in the *Object Inspector*.  These docstrings give a description of the required and optional arguments to each method, along with examples and notes where applicable.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Topology Manipulations and Operations
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
It is possible to add and remove pores and throats from the Network topology after it has been generated.  The ``trim`` command takes a list of pore or throat numbers and removes from the Network, while the 'extend' command receives a set of pore coordinates and/or throat connections and adds them to the Network:

>>> pn.trim(pores=[0,2,4])
>>> print(pn)
------------------------------------------------------------
OpenPNM.Network.Cubic: 	Cubic_2xGW2
------------------------------------------------------------
#     Properties                          Valid Values
------------------------------------------------------------
1     pore.coords                            24 / 24   
2     pore.index                             24 / 24   
3     throat.conns                           43 / 43   
------------------------------------------------------------
#     Labels                              Assigned Locations
------------------------------------------------------------
1     pore.all                            24        
2     pore.back                           9         
3     pore.bottom                         8         
4     pore.front                          6         
5     pore.internal                       24        
6     pore.left                           7         
7     pore.right                          9         
8     pore.top                            8         
9     throat.all                          43        
------------------------------------------------------------

Notice that 3 pores have indeed been removed, but also a number of throats are missing as well.  This is because throat MUST connect to a pore on both ends, so the removal of a pore necessitates the removal of all throats connected to it as well.  Throats can generally be removed without concern, however, it is very possible that isolated single pores or clusters of pores could be created that are disconnect from the main body of the network.  For instance, removing all throats connected to pore 1 will obviously lead to pore 1 being isolated from the network:

>>> Ts = pn.find_neighbor_throats(pores=1)
>>> pn.trim(throats=Ps)

The 'health' of the Network can be checked with a built-in method:

>>> pn.check_network_health()
{'duplicate_throats': [], 'isolated_pores': array([1], dtype=int64), 'disconnected_clusters': [array([ 0,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23], dtype=int64), array([1], dtype=int64)], 'bidirectional_throats': []}

The check found that pore 1 is now an 'isolated_pore'.

Extending the network can also be done.  For instance, it is possible to reconnect pore 1 to the main network:

>>> pn.extend(throat_conns=[[0,1]])
>>> pn.find_neighbor_pores(pores=1)
array([0])

This indicates that pore now has pore 0 as a connected neighbor.  A health check of the network would also pass cleanly.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Creating Custom Network Topology Generators
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For description of how to create customized subclasses, see :ref:`Customizing OpenPNM<customizing>`

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Network Topology: In Depth
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
As the name suggests, pore network modeling borrows significantly from the fields of network and graph theory.  During the development of OpenPNM, it was debated whether existing Python graph theory packages (such as `graph-tool <http://graph-tool.skewed.de/>`_ and `NetworkX <http://networkx.github.io/>`_) should be used to store the network topology.  It was decided that storage of network property data should be simply stored as 1D Numpy ndarrays.  In this form the data storage would be very transparent, since all engineers are used to working with 1D arrays (i.e. vectors), and also very efficiently since this allows a high degree of code vectorization.  Fortuitously, around the same time as this discussion, Scipy started to include the `compressed sparse graph <http://docs.scipy.org/doc/scipy/reference/sparse.csgraph.html>`_ library, which contained numerous graph theory algorithms.  The CSGraph library requires adjacency matrices which happens to be how OpenPNM stores network connections as described below.

One of the main design considerations of OpenPNM was to accommodate *all* pore networks (arbitrary dimensionality, connectivity, shape and so on).  Cubic networks are commonly used in pore network modeling, with each pore connected to 6 or 26 neighbors.  This type of network *can* be represented as cubic matrices in numerical simulations, and this has the advantage that it is easily interpreted by human users.  Representing networks this way, however, clearly lacks generality.  Networks extracted from tomographic images, or generated using random pore placements connected by Delaunay tessellations require a different approach.  OpenPNM uses network representation schemes borrowed from graph theory, such as adjacency and incidence matrices, that can be used to represent *all* network topologies.

The only topology definitions required by OpenPNM are:

1. A throat connects exactly two pores, no more and no less

2. Throats are non-directional, meaning that flow in either direction is equal (note that this restriction might be worth relaxing in a future release)

There are other general rules, but these are not essential:

3. Pores can have an arbitrary number of throats, including zero; however, pores with zero throats lead to singular matrices and other problems so should be avoided.

4. Two pores are connected by no more than one throat, unless there is some real physical reason for this.  Unintentional duplicate connections impact the rate of mass exchange between pores.  

A network has a certain number of pores, *Np*, and a certain number of throats, *Nt*.  Typically, *Nt* > *Np* since most pores have more than 1 throat.  If every pore has 1 throat (e.g. the network forms a circular chain), then *Nt* = *Np* - 1.  It can be *unofficially* stated that a network should have at least 2 pores connected by at least 1 throat (*Np* > 1 and *Nt* > 0).

-------------------------------------------------------------------------------
Storing Network Connectivity with Adjacency Matrices
-------------------------------------------------------------------------------
Network topology or connectivity is conveniently and efficiently stored as an `adjacency matrix <http://en.wikipedia.org/wiki/Adjacency_matrix>`_.  An adjacency matrix is a *Np*-by-*Np* 2D matrix.  A non-zero value at location (*i*, *j*) indicates that pores *i* and *j* are connected.  Describing the network in this fashion is one of the main features that allows OpenPNM to be agnostic to the type of network it describes.  Another important feature of the adjacency matrix is that it is highly sparse and can be stored with a variety of sparse storage schemes.  OpenPNM stores the adjacency matrix in the 'COO' or 'IJV' format, which essentially stores the coordinates (I,J) and values (V) of the nonzero elements in three separate lists.  This approach results in `throat_data` entry called *'conns'* which is and *Nt*-by-2 array that gives the ID number of the two pores that a given throat connects.  The storage scheme coincides exactly with the storage of all other throat properties.  The details of the OpenPNM implementation of adjacency matrices and other relate issues are given below.

When each pore has a unique ID number it is logical to store the network connectivity as a list of the pores to	which a given pore is connected.  Graph theoreticians have devised an elegant and powerful approach for storing this information, which OpenPNM has adopted, called adjacency matrices.  An adjacency matrix is a sparse 2D matrix of size *Np*-by-*Np*.  A value of 1 is placed at location (*i*, *j*) to indicate that pores *i* and *j* are connected.  In pore networks there is generally no difference between traversing from pore *i* to pore *j* or from pore *j* to pore *i*, so a 1 is also placed at location (*j*, *i*).  This means that determining which pores are connected directly to a given pore (say *i*) can be accomplished by finding the locations of non-zeros in row *i*.  In graph theory terminology this is deemed an *undirected* network, meaning that the *direction* of traversal is immaterial.  The adjacency matrix of an undirected network is symmetric.  Since the adjacency matrix is symmetric it is redundant to store the entire matrix when only the upper (or lower) triangular part is necessary.

Because pores are generally only connected to nearby pores, the number of throats per pore is a very small fraction of the total number of throats.  This means that there are very few non-zero elements on each row, so the adjacency matrix is highly sparse.  This fact naturally lends itself to sparse storage schemes.  OpenPNM uses uses the IJV sparse storage scheme to store the upper triangular portion of the adjacency matrix.  The *IJV* scheme is simply an *Np*-by-3 array of the (*I*, *J*) coordinates of each non-zero element in the adjacency matrix, along with the corresponding non-zero value (*V*).  (The scipy.sparse module calls this the Coordinate or COO storage scheme, but it is more widely known as IJV).  For example, to denote a value of 1 on row 3 and column 7, the *IJV* storage scheme would include an entry IJV = [3, 7, 1].  Each non-zero element in the adjacency matrix corresponds to a row to the *IJV* array.  Moreover, the number of non-zeros in the upper triangular portion of the adjacency matrix is equal to the number of throats in the network, so the dimensions of the *IJV* array is *Nt*-by-3.  This is not a coincidence; a key feature of the adjacency matrix is that each non-zero element directly corresponds to a throat.  Because throat numbers are implicitly defined by their location in an array, then the IJV sparse storage scheme automatically assigns throat ID numbers when the IJV array is generated.  For instance, when scanning the adjacency matrix from left-to-right, top-to-bottom, the first non-zero element encountered (say at location [0,5]) would be assigned throat number 0, and stored as IJV[0] = [0,5,1].

One further optimization used by OpenPNM is to drop the V from the IJV format since the non-zeros in the adjacency matrix are all 1.  This results in a *Nt*-by-2 array which is called 'throat.conns'.  Any desired throat property array can be appended as a third column to the 'throat.conns' array to fully specify the IJV format for use with the scipy.sparse or scipy.csgraph functions.  OpenPNM provides a routine for this operation (``create_adjacency_matrix``), which takes the desired throat property list to insert into *V* as an argument.

In summary, when storing network connectivity as the upper triangular portion of an adjacency in the IJV sparse storage format, the end result is an *Nt*-by-2 list describing which pores are connected by a given throat.  These connections are a fundamental property associated with each throat in the same way as throat diameter or capillary entry pressure.  This highly distilled storage format minimizes memory usage, allows for vectorization of the code, is the most efficient means of generating a sparse matrix, and corresponds perfectly with the storage of other throat properties using the ID number implicitly defined by the list element location.

**Other Sparse Storage Schemes**

The IJV storage format corresponds perfectly with the way other throat data is stored in OpenPNM, however some tasks and queries are performed more efficiently using other storage formats.  OpenPNM converts between these formats internally as needed.  For instance, most linear solvers prefer the compressed-sparse-row (CSR) scheme.  Conveniently, the IJV format used by OpenPNM is the fastest way to generate sparse matrices, so conversion, or building of each required sparse format is very efficient.  OpenPNM uses the methods provided by scipy.sparse for these conversions so they are highly optimized and based on C.  OpenPNM contains a method for constructing sparse matrices (called fill_adjacency_matrix) which accepts the storage type as an argument (i.e. 'csr', 'lil', etc).  This method can generate these other formats very quickly since they all derive from the IJV ('coo') format.  For a discussion of sparse storage schemes and the respective merits, see this `Wikipedia article <http://en.wikipedia.org/wiki/Sparse_matrix>`_.

**Incidence Matrices**

Another way to represent network connections is an incidence matrix.  This is similar to an adjacency matrix but rather than denoting which pores are connected to which, it denotes which pores are connected to which throats.  An incidence matrix is *Np*-by-*Nt* in size, with *Nt* non-zero elements.  The incidence matrix is useful for quickly querying which throats are connected to a given pore by finding the location of non-zero elements on a row.  Incidence matrices are generated as needed by OpenPNM internally for performing such queries, and the user does not usually interact with them.

