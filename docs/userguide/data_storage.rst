.. _data_storage:

###############################################################################
Network Architecture and Data Storage
###############################################################################

OpenPNM utilizes the object oriented capacities of Python by defining a network as an object.  A network object contains both the data that describes the network along with the tools, functions, and methods needed to access this data in ways applicable to the pore network modeling paradigm.  One key feature of this object is that it is completely agnostic about the type of network it describes; a random, cubic or another network topology is stored in exactly the same manner.
The advantages of this are numerous:

1. The OpenPNM framework can be applied to any network or situation
2. All operations performed on the network can be fully generic
3. Only one version of an algorithm or method needs to be developed then can be applied universally
4. The commonality between network types becomes apparent

As the name suggests, pore network modeling borrows significantly from the fields of network and graph theory.  During the development of OpenPNM, it was debated whether existing Python graph theory packages (such as `graph-tool <http://graph-tool.skewed.de/>`_ and `NetworkX <http://networkx.github.io/>`_) should be used to store the network topology.  It was decided that storage of network property data could be handled very efficiently and transparently using simple 1D arrays, which allows for a high degree of code vectorization.  Fortuitously, around the same time as this discussion, Scipy started to include the 'compressed sparse graph' library, which contained numerous graph theory algorithms.  The CSGraph library requires adjacency matrices which happens to be how OpenPNM stores network connections as described below.

===============================================================================
Network Architecture
===============================================================================

One of the main design considerations of OpenPNM was to accommodate *all* pore networks (arbitrary dimensionality, connectivity, shape and so on).  Cubic networks are commonly used in pore network modeling, with each pore connected to 6 or 26 neighbors.  This type of network *can* be represented as cubic matrices in numerical simulations, and this has the advantage that it is easily interpreted by human users.  Representing networks this way, however, clearly lacks generality.
Networks extracted from tomographic images, or generated using random pore placements connected by Delaunay tessellations require a different approach.  OpenPNM uses network representation schemes borrowed from graph theory, such as adjacency and incidence matrices, that can be used to represent *all* network topologies.

The basic definitions used by OpenPNM are:

1. Pores can have an arbitrary number of throats greater than zero

2. A throat connects exactly two pores, no more and no less

3. Two pores are connected by no more than one throat

A network has a certain number of pores, *Np*, and a certain number of throats, *Nt*.  Typically, *Nt* > *Np* since most pores have more than 1 throat.  If every pore has 1 throat (e.g. the network forms a circular chain), then *Nt* = *Np* - 1.
It can be *unofficially* stated that a network should have at least 2 pores connected by at least 1 throat (*Np* > 1 and *Nt* > 0).

-------------------------------------------------------------------------------
Storing Network Connectivity with Adjacency Matrices
-------------------------------------------------------------------------------

Network topology or connectivity is conveniently and efficiently stored as an `adjacency matrix <http://en.wikipedia.org/wiki/Adjacency_matrix>`_.  An adjacency matrix is a *Np*-by-*Np* 2D matrix.  A non-zero value at location (*i*, *j*) indicates that pores *i* and *j* are connected.  Describing the network in this fashion is one of the main features that allows OpenPNM to be agnostic to the type of network it describes.  Another important feature of the adjacency matrix is that it is highly sparse and can be stored with a variety of sparse storage schemes.  OpenPNM stores the adjacency matrix in the 'COO' or 'IJV' format, which essential stores the coordinates (I,J) and values (V) of the nonzero elements.  Without delving into the details, this approach results in `throat_properties` entry called *'connections'* which is and *Nt*-by-2 array that gives the ID number of the two pores that a given throat connects.  The storage scheme coincides exactly with the storage of all other throat properties.  The details of the OpenPNM implementation of adjacency matrices and other relate issues are given below for the interested reader.

.. Note:: In Depth: Adjacency and Incidence Matrices

	*Adjacency Matrices*

	When each pore has a unique ID number it is logical to store the network connectivity as a list of the pores to
	which a given pore is connected.  Graph theoreticians have devised an elegant and powerful approach for storing this information, which OpenPNM has adopted, called adjacency matrices.  An adjacency matrix is a sparse 2D matrix of size *Np*-by-*Np*.  A value of 1 is placed at location (*i*, *j*) to indicate that pores *i* and *j* are connected.  In pore networks there is generally no difference between traversing from pore *i* to pore *j* or from pore *j* to pore *i*, so a 1 is also placed at location (*j*, *i*).  This means that determining which pores are connected directly to a given pore (say *i*) can be accomplished by finding the locations of non-zeros in row *i*.  In graph theory terminology this is deemed an *undirected* network, meaning that the *direction* of traversal is immaterial.  The adjacency matrix of an undirected network is symmetric.  Since the adjacency matrix is symmetric it is redundant to store the entire matrix when only the upper (or lower) triangular part is necessary.

	Because pores are generally only connected to nearby pores, the number of throats per pore is a very small fraction of the total number of throats.  This means that there are very few non-zero elements on each row, so the adjacency matrix is highly sparse.  This fact naturally lends itself to sparse storage schemes.  OpenPNM uses uses the IJV sparse storage scheme to store the upper triangular portion of the adjacency matrix.  The *IJV* scheme is simply an *Np*-by-3 array of the (*I*, *J*) coordinates of each non-zero element in the adjacency matrix, along with the corresponding non-zero value (*V*).  (The scipy.sparse module calls this the Coordinate or COO storage scheme, but it is more widely known as IJV).  For example, to denote a value of 1 on row 3 and column 7, the *IJV* storage scheme would include an entry IJV = [3, 7, 1].  Each non-zero element in the adjacency matrix corresponds to a row to the *IJV* array.  Moreover, the number of non-zeros in the upper triangular portion of the adjacency matrix is equal to the number of throats in the network, so the dimensions of the *IJV* array is *Nt*-by-3.  This is not a coincidence; a key feature of the adjacency matrix is that each non-zero element directly corresponds to a throat.  Because throat numbers are implicitly defined by their location in an array, then the IJV sparse storage scheme automatically assigns throat ID numbers when the IJV array is generated.  For instance, when scanning the adjacency matrix from left-to-right, top-to-bottom, the first non-zero element encountered (say at location [0,5]) would be assigned throat number 0, and stored as IJV[0] = [0,5,1].

	One further optimization used by OpenPNM is to drop the V from the IJV format since the non-zeros in the adjacency matrix are all 1.  This results in a *Nt*-by-2 array which is called *connections*.  Any desired throat property array can be appended as a third column to the *connections* array to fully specify the IJV format for use with the scipy.sparse or scipy.csgraph functions.  OpenPNM provides a routine for this operation (``'fill_adjacency_matrix'``), which takes the desired throat property list to insert into *V* as an argument.

	In summary, when storing network connectivity as the upper triangular portion of an adjacency in the IJV sparse storage format, the end result is an *Nt*-by-2 list describing which pores are connected by a given throat.  These connections are a fundamental property associated with each throat in the same way as throat diameter or capillary entry pressure.  This highly distilled storage format minimized memory usage, allows for vectorization of the code, is the most efficient means of generating a sparse matrix, and corresponds perfectly with the storage of other throat properties using the ID number implicitly defined by the list element location.

	*Other Sparse Storage Schemes*

	The IJV storage format corresponds perfectly with the way other throat data is stored in OpenPNM, however some tasks and queries are performed more efficiently using other storage formats.  OpenPNM converts between these formats internally as needed.  For instance, most linear solvers prefer the compressed-sparse-row (CSR) scheme.  Conveniently, the IJV format used by OpenPNM is the fastest way to generate sparse matrices, so conversion, or building of each required sparse format is very efficient.  OpenPNM uses the methods provided by scipy.sparse for these conversions so they are highly optimized and based on C.  OpenPNM contains a method for constructing sparse matrices (called fill_adjacency_matrix) which accepts the storage type as an argument (i.e. 'csr', 'lil', etc).  This method can generate these other formats very quickly since they all derive from the IJV ('coo') format.  For a discussion of sparse storage schemes and the respective merits, see this `Wikipedia article <http://en.wikipedia.org/wiki/Sparse_matrix>`_.

	*Incidence Matrices*

	Another way to represent network connections is an incidence matrix.  This is similar to an adjacency matrix but rather than denoting which pores are connected to which, it denotes which pores are connected to which throats.  An incidence matrix is *Np*-by-*Nt* in size, with *Nt* non-zero elements.  The incidence matrix is useful for quickly querying which throats are connected to a given pore by finding the location of non-zero elements on a row.  Incidence matrices are generated as needed by OpenPNM internally for performing such queries, and the user does not usually interact with them.

===============================================================================
Data Storage
===============================================================================

OpenPNM stores all data in 1D arrays or lists.  This format is well suited for vectorized calculations which are essential for fast and efficient computations (see Note below).  Storing data as 1D lists also allows for a topologically agnostic network framework, since cubic and random networks are all stored in the same list format.  As discussed above, the connectivity in the network is tracked using adjacency matrices.  Storage of all data in 1D lists means that each pore (or throat) is implicitly assigned an ID number, which corresponds to it's location in the list.  Specifically, if list A contains pore diameter and list B contains pore volume, then `A[6]` is the diameter of pore `6` and and `B[6]` contains it's volume.  

.. Note:: Numpy ND-arrays
   
   OpenPNM stores all pore and throat properties as Numpy ndarrays.  ndarrays are a numerical data type provided by the Numpy package (which is embedded in the Scipy package) that allow for the type of numerical manipulations that scientists and engineers expect, such as vectorization, slicing, boolean indexing and so on.

Another important aspect of the data storage scheme is that pore and throat data are stored separately.  This is to prevent properties with the same name from colliding (such as volume).  OpenPNM uses the Python dictionary data-type to store each property by name, either in the pore_data or throat_data dictionary.  For instance, pore volumes are stored as pore_data['volume'], while throat volumes are stored as throat_data['volume'].  This approach ensures that all data stored in the same dictionary are of the same length (*Nt* or *Np*).  


-------------------------------------------------------------------------------
Pore and Throat *Data* and *Info*
-------------------------------------------------------------------------------
OpenPNM stores two types of information about pores and throats: 'data' and 'info'.  Data includes the physical quantities associated with a pore or throat such as the geometrical (e.g. diameter), structural (e.g. coordinates) and thermophysical (e.g. capillary entry pressure) aspects of the network.  Info is basically pore or throat labels, such as which subdomain a pore belongs to, whether a pore is an internal or boundary pore, and so on. 

`Data` is 




-------------------------------------------------------------------------------
Data and Info: Setter and Getter Methods
-------------------------------------------------------------------------------




















