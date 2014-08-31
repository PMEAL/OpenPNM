.. _data_storage:

###############################################################################
Network Architecture and Data Storage
###############################################################################
As the name suggests, pore network modeling borrows significantly from the fields of network and graph theory.  During the development of OpenPNM, it was debated whether existing Python graph theory packages (such as `graph-tool <http://graph-tool.skewed.de/>`_ and `NetworkX <http://networkx.github.io/>`_) should be used to store the network topology.  It was decided that storage of network property data should be simply stored as 1D Numpy ndarrays.  In this form the data storage would be very transparent, since all engineers are used to working with 1D arrays (i.e. vectors), and also very efficiently since this allows a high degree of code vectorization.  Fortuitously, around the same time as this discussion, Scipy started to include the `compressed sparse graph <http://docs.scipy.org/doc/scipy/reference/sparse.csgraph.html>`_ library, which contained numerous graph theory algorithms.  The CSGraph library requires adjacency matrices which happens to be how OpenPNM stores network connections as described below.

===============================================================================
Data Storage
===============================================================================
Before proceeding with an overview of the framework, a few words must be said about how OpenPNM stores data, the type of data it stores, and some basic rules it applies to data formats and names. 

1. All pore and throat data are stored as lists (actually Numpy arrays).  If a network has 100 pores then all pore data (i.e. diameter) is stored in a list 100 elements long.  Each pore and throat is associated with a number, which corresponds to the location in the pore or throat lists where its information is stored (i.e. data for pore 10 is stored in element 10).  This list format allows any topology to be described by the framework, and also enables easy vectorization of code for performance gains.  A detailed description of this storage scheme is given :ref:`here <data_storage>`. 

2. OpenPNM differentiates between two types of data for pores and throats.  The physical details about pores and throats is referred as *properties*, and this includes information such as pore volume and throat length.  The second type of information is referred to as *labels*.  Labels were conceived as a means to dynamically create groups of pores and throats so they could be quickly accessed by the user.  For instance, in a cubic network it is helpful to know which pores are on the 'top' surface.  This label is automatically added by the topology generator, so a list of all pores on the 'top' can be retrieved by simply querying which pores possess the label 'top'. 

3. OpenPNM forces all data and label *names* to be either of the form 'pore.name' or 'throat.name'.  All 'pore.name' lists are Np long, and all 'throat.name' lists are Nt long, where Np and Nt are the number of pores and throats.  Forcing all arrays to be of the same length ensures that vectorized code operations will always work, and enforcing this naming allows the framework to make each array to correct length.  It is also very convenient to see what type of data is stored in a list just by glancing at its name.  

===============================================================================
Network Architecture
===============================================================================
One of the main design considerations of OpenPNM was to accommodate *all* pore networks (arbitrary dimensionality, connectivity, shape and so on).  Cubic networks are commonly used in pore network modeling, with each pore connected to 6 or 26 neighbors.  This type of network *can* be represented as cubic matrices in numerical simulations, and this has the advantage that it is easily interpreted by human users.  Representing networks this way, however, clearly lacks generality.  Networks extracted from tomographic images, or generated using random pore placements connected by Delaunay tessellations require a different approach.  OpenPNM uses network representation schemes borrowed from graph theory, such as adjacency and incidence matrices, that can be used to represent *all* network topologies.

The basic definitions used by OpenPNM are:

1. Pores can have an arbitrary number of throats greater than zero

2. A throat connects exactly two pores, no more and no less

3. Two pores are connected by no more than one throat

A network has a certain number of pores, *Np*, and a certain number of throats, *Nt*.  Typically, *Nt* > *Np* since most pores have more than 1 throat.  If every pore has 1 throat (e.g. the network forms a circular chain), then *Nt* = *Np* - 1.  It can be *unofficially* stated that a network should have at least 2 pores connected by at least 1 throat (*Np* > 1 and *Nt* > 0).

-------------------------------------------------------------------------------
Storing Network Connectivity with Adjacency Matrices
-------------------------------------------------------------------------------

Network topology or connectivity is conveniently and efficiently stored as an `adjacency matrix <http://en.wikipedia.org/wiki/Adjacency_matrix>`_.  An adjacency matrix is a *Np*-by-*Np* 2D matrix.  A non-zero value at location (*i*, *j*) indicates that pores *i* and *j* are connected.  Describing the network in this fashion is one of the main features that allows OpenPNM to be agnostic to the type of network it describes.  Another important feature of the adjacency matrix is that it is highly sparse and can be stored with a variety of sparse storage schemes.  OpenPNM stores the adjacency matrix in the 'COO' or 'IJV' format, which essential stores the coordinates (I,J) and values (V) of the nonzero elements.  Without delving into the details, this approach results in `throat_data` entry called *'conns'* which is and *Nt*-by-2 array that gives the ID number of the two pores that a given throat connects.  The storage scheme coincides exactly with the storage of all other throat properties.  The details of the OpenPNM implementation of adjacency matrices and other relate issues are given below for the interested reader.

*Adjacency Matrices*

When each pore has a unique ID number it is logical to store the network connectivity as a list of the pores to	which a given pore is connected.  Graph theoreticians have devised an elegant and powerful approach for storing this information, which OpenPNM has adopted, called adjacency matrices.  An adjacency matrix is a sparse 2D matrix of size *Np*-by-*Np*.  A value of 1 is placed at location (*i*, *j*) to indicate that pores *i* and *j* are connected.  In pore networks there is generally no difference between traversing from pore *i* to pore *j* or from pore *j* to pore *i*, so a 1 is also placed at location (*j*, *i*).  This means that determining which pores are connected directly to a given pore (say *i*) can be accomplished by finding the locations of non-zeros in row *i*.  In graph theory terminology this is deemed an *undirected* network, meaning that the *direction* of traversal is immaterial.  The adjacency matrix of an undirected network is symmetric.  Since the adjacency matrix is symmetric it is redundant to store the entire matrix when only the upper (or lower) triangular part is necessary.

Because pores are generally only connected to nearby pores, the number of throats per pore is a very small fraction of the total number of throats.  This means that there are very few non-zero elements on each row, so the adjacency matrix is highly sparse.  This fact naturally lends itself to sparse storage schemes.  OpenPNM uses uses the IJV sparse storage scheme to store the upper triangular portion of the adjacency matrix.  The *IJV* scheme is simply an *Np*-by-3 array of the (*I*, *J*) coordinates of each non-zero element in the adjacency matrix, along with the corresponding non-zero value (*V*).  (The scipy.sparse module calls this the Coordinate or COO storage scheme, but it is more widely known as IJV).  For example, to denote a value of 1 on row 3 and column 7, the *IJV* storage scheme would include an entry IJV = [3, 7, 1].  Each non-zero element in the adjacency matrix corresponds to a row to the *IJV* array.  Moreover, the number of non-zeros in the upper triangular portion of the adjacency matrix is equal to the number of throats in the network, so the dimensions of the *IJV* array is *Nt*-by-3.  This is not a coincidence; a key feature of the adjacency matrix is that each non-zero element directly corresponds to a throat.  Because throat numbers are implicitly defined by their location in an array, then the IJV sparse storage scheme automatically assigns throat ID numbers when the IJV array is generated.  For instance, when scanning the adjacency matrix from left-to-right, top-to-bottom, the first non-zero element encountered (say at location [0,5]) would be assigned throat number 0, and stored as IJV[0] = [0,5,1].

One further optimization used by OpenPNM is to drop the V from the IJV format since the non-zeros in the adjacency matrix are all 1.  This results in a *Nt*-by-2 array which is called *connections*.  Any desired throat property array can be appended as a third column to the *connections* array to fully specify the IJV format for use with the scipy.sparse or scipy.csgraph functions.  OpenPNM provides a routine for this operation (``'fill_adjacency_matrix'``), which takes the desired throat property list to insert into *V* as an argument.

In summary, when storing network connectivity as the upper triangular portion of an adjacency in the IJV sparse storage format, the end result is an *Nt*-by-2 list describing which pores are connected by a given throat.  These connections are a fundamental property associated with each throat in the same way as throat diameter or capillary entry pressure.  This highly distilled storage format minimized memory usage, allows for vectorization of the code, is the most efficient means of generating a sparse matrix, and corresponds perfectly with the storage of other throat properties using the ID number implicitly defined by the list element location.

*Other Sparse Storage Schemes*

The IJV storage format corresponds perfectly with the way other throat data is stored in OpenPNM, however some tasks and queries are performed more efficiently using other storage formats.  OpenPNM converts between these formats internally as needed.  For instance, most linear solvers prefer the compressed-sparse-row (CSR) scheme.  Conveniently, the IJV format used by OpenPNM is the fastest way to generate sparse matrices, so conversion, or building of each required sparse format is very efficient.  OpenPNM uses the methods provided by scipy.sparse for these conversions so they are highly optimized and based on C.  OpenPNM contains a method for constructing sparse matrices (called fill_adjacency_matrix) which accepts the storage type as an argument (i.e. 'csr', 'lil', etc).  This method can generate these other formats very quickly since they all derive from the IJV ('coo') format.  For a discussion of sparse storage schemes and the respective merits, see this `Wikipedia article <http://en.wikipedia.org/wiki/Sparse_matrix>`_.

*Incidence Matrices*

Another way to represent network connections is an incidence matrix.  This is similar to an adjacency matrix but rather than denoting which pores are connected to which, it denotes which pores are connected to which throats.  An incidence matrix is *Np*-by-*Nt* in size, with *Nt* non-zero elements.  The incidence matrix is useful for quickly querying which throats are connected to a given pore by finding the location of non-zero elements on a row.  Incidence matrices are generated as needed by OpenPNM internally for performing such queries, and the user does not usually interact with them.











