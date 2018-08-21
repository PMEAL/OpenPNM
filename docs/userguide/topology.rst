.. _topology:

================================================================================
Representing Topology
================================================================================

.. contents:: Page Contents
    :depth: 3

--------------------------------------------------------------------------------
Storage of Topological Connections
--------------------------------------------------------------------------------
As the name suggests, pore network modeling borrows significantly from the fields of network and graph theory.  During the development of OpenPNM, it was debated whether existing Python graph theory packages (such as `graph-tool <http://graph-tool.skewed.de/>`_ and `NetworkX <http://networkx.github.io/>`_) should be used to store the network topology.  It was decided that network property data should be simply stored as `Numpy ND-arrays <http://www.numpy.org/>`_).  This format makes the data storage very transparent and familiar since all engineers are used to working with arrays (i.e. vectors), and also very efficiently since this allows code vectorization.  Fortuitously, around the same time as this discussion, Scipy introduced the `compressed sparse graph library <http://docs.scipy.org/doc/scipy/reference/sparse.csgraph.html>`_, which contains numerous graph theory algorithms that take Numpy arrays as arguments.  Therefore, OpenPNM's topology model is implemented using Numpy arrays, which is described in detail below:

The only topology definitions required by OpenPNM are:

1. A throat connects exactly two pores, no more and no less

2. Throats are non-directional, meaning that flow in either direction is equal

Other general, but non-essential rules are:

3. Pores can have an arbitrary number of throats, including zero; however, pores with zero throats lead to singular matrices and other problems so should be avoided.

4. Two pores are generally connected by no more than one throat.  It is technically possible in OpenPNM to have multiple throats between a pair of pores, but it is not rigorosly supported so unintended results may arise.

-------------------------------------------------------------------------------
Sparse Adjacency Matrices
-------------------------------------------------------------------------------

In OpenPNM network topology (or connectivity) is stored as an `adjacency matrix <http://en.wikipedia.org/wiki/Adjacency_matrix>`_.  An adjacency matrix is a *Np*-by-*Np* 2D matrix.  A non-zero value at location (*i*, *j*) indicates that pores *i* and *j* are connected.  Describing the network in this general fashion allows OpenPNM to be agnostic to the type of network it describes.  Another important feature of the adjacency matrix is that it is highly sparse and can be stored with a variety of sparse storage schemes.  OpenPNM stores the adjacency matrix in the 'COO' or 'IJV' format, which essentially stores the coordinates (I,J) and values (V) of the nonzero elements in three separate lists.  This approach results in a property called ``'throat.conns'``; it is an *Nt*-by-2 array that gives the index of the two pores on either end of a given throat.  The representation of an arbitrary network is shown in the following figure. It has 5 pores and 7 throats, and the ``'throat.conns'`` array contains the (I,J,V) information to describes the adjacency matrix.

.. image:: http://i.imgur.com/rMpezCc.png
    :width: 500 px
    :align: center

.. note:: Additional Thoughts on Sparse Storage

    * In pore networks there is generally no difference between traversing from pore *i* to pore *j* or from pore *j* to pore *i*, so a 1 is also found at location (*j*, *i*) and the matrix is symmetrical.

    * Since the adjacency matrix is symmetric, it is redundant to store the entire matrix when only the upper triangular part is necessary.  The ``'throat.conns'`` array only stores the upper triangular information, and *i* is always less than *j*.

    * Although this storage scheme is widely known as *IJV*, the ``scipy.sparse`` module calls this the Coordinate or *COO* storage scheme.

    * Some tasks are best performed on other types of storages scheme, such as *CSR* or *LIL*.  OpenPNM converts between these internally as necessary, but users can generate a desired format using the ``create_adjacency_matrix`` method which accepts the storage type as an argument (i.e. ``'csr'``, ``'lil'``, etc).  For a discussion of sparse storage schemes and the respective merits, see this `Wikipedia article <http://en.wikipedia.org/wiki/Sparse_matrix>`_.

--------------------------------------------------------------------------------
Performing Network Queries
--------------------------------------------------------------------------------

Querying and inspecting the pores and throats in the **Network** is an important tool for working with networks. The various functions that are included on the **GenericNetwork** class will be demonstrated below on the following cubic network:

.. code-block:: python

    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[10, 10, 10])

................................................................................
Finding Neighboring Pores and Throats
................................................................................

Given a pore *i*, it possible to find which pores (or throats) are directly connected to it:

.. code-block:: Python

    >>> Ps = pn.find_neighbor_pores(pores=1)
    >>> print(Ps)
    [  0   2  11 101]
    >>> Ps = pn.find_neighbor_throats(pores=1)
    >>> print(Ps)
    [   0    1  901 1801]

The above queries can be more complex if a list of pores is sent, and the ```mode``` argument is specified.  This is useful for finding neighbors surrounding a set of pores such as the fringes around an invading fluid cluster, or all throats within a cluster:

.. code-block:: python

    >>> Ps = pn.find_neighbor_pores(pores=[2, 3, 4], mode='or')  # 'union' is default
    >>> print(Ps)
    [  1   5  12  13  14 102 103 104]
    >>> Ts = pn.find_neighbor_throats(pores=[2, 3, 4], mode='xnor')
    >>> print(Ts)
    [2 3]
    >>> Ts = pn.find_neighbor_throats(pores=[2, 3, 4], mode='exclusive_or')
    >>> print(Ts)
    [   1    4  902  903  904 1802 1803 1804]

The ```mode``` argument limits the returned results using *set-theory* type logic.  Consider the following two queries:

.. code-block:: python

    >>> Ts = pn.find_neighbor_throats(pores=2)
    >>> print(Ts)
    [   1    2  902 1802]
    >>> Ts = pn.find_neighbor_throats(pores=3)
    >>> print(Ts)
    [   2    3  903 1803]

The *or* is a single set of unique values obtained by combining the two sets, while the *intersection* of these two sets includes only the values present in both (i.e. *2*)  The *difference* of these sets is all the values except those found common to both initial sets.  It's possible to specify as many pores as desired, and the *set-logic* is bit less obvious.  More generally:

* ``'or'`` returns a list of unique locations neighboring any input pores
* ``'xor'`` returns a list of locations that are only neighbors to one of the input pores
* ``'xnor'`` returns a list of locations that are neighbors to at least two inputs pores

In addition to these neighbor lookups, the GenericNetwork class also offers several other methods that complete the suite of lookup tools:  ``find_connected_pores``, ``find_connecting_throats`` and ``find_nearby_pores``.  There are also many more tools related to Network queries and manipulations in the :ref:`topotools_index` module.
