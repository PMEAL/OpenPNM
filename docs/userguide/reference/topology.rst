.. _topology:

###############################################################################
Representing Topology in OpenPNM
###############################################################################

===============================================================================
Storage of Topological Connections
===============================================================================
As the name suggests, pore network modeling borrows significantly from the fields of network and graph theory.  During the development of OpenPNM, it was debated whether existing Python graph theory packages (such as `graph-tool <http://graph-tool.skewed.de/>`_ and `NetworkX <http://networkx.github.io/>`_) should be used to store the network topology.  It was decided that storage of network property data should be simply stored as Numpy ndarrays (see `Numpy <http://www.numpy.org/>`_).  In this form the data storage would be very transparent, since all engineers are used to working with arrays (i.e. vectors), and also very efficiently since this allows code vectorization.  Fortuitously, around the same time as this discussion, Scipy introduced the `compressed sparse graph library <http://docs.scipy.org/doc/scipy/reference/sparse.csgraph.html>`_, which contains numerous graph theory algorithms that take Numpy arrays as arguments.

The only topology definitions required by OpenPNM are:

1. A throat connects exactly two pores, no more and no less

2. Throats are non-directional, meaning that flow in either direction is equal (note that this restriction might be worth relaxing in a future release)

Other general, but non-essential rules are:

3. Pores can have an arbitrary number of throats, including zero; however, pores with zero throats lead to singular matrices and other problems so should be avoided.

4. Two pores are connected by no more than one throat, unless there is some real physical reason for this.  Unintentional duplicate connections impact the rate of mass exchange between pores.

The **GenericNetwork** class has a ``check_network_health`` method that scans the network for the above criteria as well as a few others and returns a **HealthDict** which lists if any problems were founds, and where.

-------------------------------------------------------------------------------
Sparse Adjacency Matrices
-------------------------------------------------------------------------------

In OpenPNM network topology (or connectivity) is stored as an `adjacency matrix <http://en.wikipedia.org/wiki/Adjacency_matrix>`_.  An adjacency matrix is a *Np*-by-*Np* 2D matrix.  A non-zero value at location (*i*, *j*) indicates that pores *i* and *j* are connected.  Describing the network in this general fashion allows OpenPNM to be agnostic to the type of network it describes.  Another important feature of the adjacency matrix is that it is highly sparse and can be stored with a variety of sparse storage schemes.  OpenPNM stores the adjacency matrix in the 'COO' or 'IJV' format, which essentially stores the coordinates (I,J) and values (V) of the nonzero elements in three separate lists.  This approach results in a property which OpenPNM calls ``'throat.conns'``; it is an *Nt*-by-2 array that gives the ID number of the two pores on either end of a given throat.  The representation of an arbitrary network is shown in following figure. It has 5 pores and 7 throats, and the ``'throat.conns'`` array contains the (I,J,V) information to describes the adjacency matrix.

.. image:: http://i.imgur.com/rMpezCc.png

-------------------------------------------------------------------------------
Additional Thoughts on Topology Storage
-------------------------------------------------------------------------------

* In pore networks there is generally no difference between traversing from pore *i* to pore *j* or from pore *j* to pore *i*, so a 1 is also found at location (*j*, *i*) and the matrix is symmetrical.

* A symmetrical matrix means that determining which pores are connected directly to a given pore (say *i*) can be accomplished by finding the locations of non-zeros in row *i*.

* Since the adjacency matrix is symmetric, it is redundant to store the entire matrix when only the upper triangular part is necessary.  The ``'throat.conns'`` array only stores the upper triangular information, and *i* is always less than *j*.

* Although this storage scheme is widely known as *IJV*, the ``scipy.sparse`` module calls this the Coordinate or *COO* storage scheme.

* Some tasks are best performed on other types of storages scheme, such as *CSR* or *LIL*.  OpenPNM converts between these internally as necessary, but uses can generate a desired format using the ``create_adjacency_matrix`` method which accepts the storage type as an argument (i.e. ``'csr'``, ``'lil'``, etc).  For a discussion of sparse storage schemes and the respective merits, see this `Wikipedia article <http://en.wikipedia.org/wiki/Sparse_matrix>`_.

===============================================================================
Performing Network Queries
===============================================================================



===============================================================================
Manipulating and Altering Topology
===============================================================================
