===============================================================================
Network Architecture and Data Storage Formats
===============================================================================
OpenPNM utilizes the object oriented capacities of Python.  The code is built upon the idea that a network is an object.  A network object contains both the data that describes the network properties, and the tools, functions, or methods needed to access this data in ways applicable to the pore network modeling paradigm.  One key feature of this object is that it is completely agnostic about the type of network it describes; a random, cubic or another netowrk topology is stored in exactly the same manner.  The most important repercussion of this choice is the fact that all physical algorithms (such as diffusion or drainage) operating on the network can be fully generic, and the all methods that read and write network property data can be fully generic as well.  

As the name suggests, pore network modeling borrows signifcantly from the fields of network and graph theory.  During the development of OpenPNM, it was debated whether existing Python graph theory packages (such as `graph-tool <http://graph-tool.skewed.de/>`_ and `NetworkX <http://networkx.github.io/>`_) should be used to store the network topology.  It was decided that storage of network property data could be handled very efficiently using 1D arrays, which allowed for a high degree of code vectorization.  Fortuitously, around the same time as this dicussion, Scipy started to include the 'compressed sparse graph' library, which contained numerous graph theory algorithms.  The CSGraph libary requires adjacency matrices in a compressed sparse storage scheme, which happens to be how OpenPNM stores network connections as descrbed below.

-------------------------------------------------------------------------------
Network Architecture
-------------------------------------------------------------------------------

*******************************************************************************
Pore and Throat Numbering
*******************************************************************************

One of the main design principle of OpenPNM was to accomodate pore networks of arbitrary dimensionality and shape.  Cubic networks are commonly used in pore network modeling, with each pore connected to 6 or 26 neighbors.  The type of network *can* be represented as cubic matrices in numerical simulations, and this has the advantage that it is easily interpreted by human users.  Representing networks this way, however, clearly lacks generality.  Networks extracted from tomographic images, or generated using random pore placements connected by Delaunnay tesselations require a different approach.  OpenPNM uses network representation schemes borrowed from graph theory, such as adjacency and incidence matrices, that can be used to represent *all* network topologies. 

The basic definitions used by OpenPNM are:

1. Pores can have an arbitrary number of throats, including zero

2. A throat connects exactly two pores, no more and no less

3. Two pores are connected by no more than one throat

A network has a certain number of pores, *Np*, and a certain number of throats, *Nt*.  Typically, *Nt* > *Np* since most pores have more than 1 throat.  If every pore has 1 throat (i.e. the network forms a linear chain), then *Nt* = *Np* - 1.  Of course, Nt can be zero but this would not be a useful network.  It can be unofficially stated that a network should have at least 2 pores connected by at least 1 throat (*Np* > 1 and *Nt * > 0).  

Each pore and throat in the network has a unique ID number.  (Python uses 0-based array indexing so the ID numbers start at 0).  In OpenPNM the ID number is *implied* by array element number, meaning that any information stored in element 0 of an array implicity applies to pore (or throat) 0.  Or in other words, finding information about pore (or throat) 0 is accomplished by looking into element 0 of an array.  There is no correspondence between pore number and throat number, meaning that throat *i* may or may not connected with pore *i*.

*******************************************************************************
Adjacency Matrices
*******************************************************************************

When each pore has a unique ID number it is logical to store the network connectivity as a list pore ID numbers to which a given pore is connected.  Although it is possible to store these *lists* as a list-of-*lists*, graph theoreticians have devised an more elegant and powerful approach, which OpenPNM has adopted, called adjacency matrices.  An adjacency matrix is a 2D matrix of size *Np*-by-*Np*.  A value of 1 is placed at location (*i*, *j*) to indicate that pores *i* and *j* are connected.  In pore networks there is generally no difference between traversing from pore *i* to pore *j* or from pore *j* to pore *i*, so a 1 is also placed at location (*j*, *i*).  This means that determining which pores are connected directly to a given pore (say *i*) can be accomplished by finding the locations of non-zeros in row *i*.  In graph theory terminology this is deemed an *undirected* network, meaning that the *direction* of traversal is immaterial.  The adjacency matrix of an undirected network is therefore symetric.  

Of course, since the adjacency matrix is symetric it is redundant to store the entire matrix when only the upper triangular part is necessary.  Furthermore, because pores are generally only connected to nearby pores, the number of throats per pore is a very small faction of the total number of throats.  This means that there are very few non-zero elements on each row, so the adjacency matrix is highly sparse.  This fact naturally lends itself to sparse storage schemes.  OpenPNM uses uses teh IJV sparse storage scheme to store the upper triangular poriton of the adjacency matrix.  The *IJV* scheme is simply an *N*-by-3 array of the (*I*, *J*) coordinates of each nonzero element in the adjacency matrix, along with the corresponding non-zero value (*V*).  (The scipy.sparse module calls this the Coordinate or COO storage scheme, but it is more widely known as IJV).  For example, to denote a value of 1 on row 3 and column 7, the *IJV* storage scheme would include an entry IJV = [3, 7, 1].  Each nonzero element in the adjacency matrix corresponds to a row to the *IJV* array.  Moreover, the number of non-zeros in the upper triangular portion of the adjacency matrix is equal to the number of throats in the network, so the dimensions of the *IJV* array is *Nt*-by-3.  This is not a coincidence; a key feature of the adjacency matrix is that each non-zero element direclty corresponds to a throat.  Because throat numbers are implicitly defined by their location in an array, then the IJV sparse storage scheme automatically assigns throat ID numbers when the IJV array is generated.  For instance, when scanning the adjacency matrix from left-to-right, top-to-bottom, the first non-zero element encountered (say at location [0,5]) would be assigned throat number 0, and stored as IJV[0] = [0,5,1].  

One further optimization used by OpenPNM is to drop the V from the IJV format since the non-zeros in the adjacency matrix are all 1.  This results in a *Nt*-by-2 array which is called *connections*.  Any desired throat property array can be appended as a third column to the *connections* array to fully specify the IJV format for use with the scipy.sparse or scipy.csgraph functions.  OpenPNM provides a routine for this operation ('fill_adjacency_matrix'), which takes the desired throat property list to insert into *V* as an argument.  

In summary, when storing network connectivity as the upper triangular portion of an adjacency in the IJV sparse storage format, the end result is an *Nt*-by-2 list describing which pores are connected by a given throat.  These connections are a fundamental property associated with each throat in the same way as throat diameter or capillary entry pressure.  This highly distilled storage format minimized memory usage, allows for vectorization of the code, is the most efficient means of generating a sparse matrix, and corresponds perfectly with the storage of other throat properties using the ID number implicitly defined by the list element location. 

*******************************************************************************
Other Sparse Storage Schemes
*******************************************************************************
The IJV storage format corresponds perfectly with the way other throat data is stored in OpenPNM, however some tasks and queries are performed more efficiently using other storage formats.  OpenPNM converts between these formats internally as needed.  

*******************************************************************************
Incidence Matrices
*******************************************************************************
Another way to represent network connections is an incidence matrix.  This is similar to an adjacency matrix but rather than denoting which pores are connected to which, it denotes which pores are connected to which throats.  An incidence matrix is *Np*-by-*Nt* is size, with *Nt* nonzero elements.  The incidence matrix is useful for quickly querying which throats are connected to a given pore by finding the location of non-zero elements on a row.  Incidence matrices are generated as needed by OpenPNM internally for performing such queries, and the user does not usually interact with them.  

-------------------------------------------------------------------------------
Network Data Storage
-------------------------------------------------------------------------------
OpenPNM stores all pore and throat properties as Numpy ndarrays (a Numpy numerical data type).  Pore properties are stored as arrays of size Np x n, where Np is the number of pores in the network and n is almost always 1, (e.g. pore volume is stored as an *Np*-by-1 array), with a few expections (e.g. spatial coordinates are stored as *Np*-by-3 for 3-dimensional space).  Throat properties are almost always stored as Nt x m arrays where Nt is the number of throats in the network.  Again, m is almost always 1 with a notable exception being the connections property that is discussed in detail below. 

**For both pores and throats, the property for pore or throat i is stored in element i of the corresponding property array.**

Assuming a pore network called 'pn', the pore properties are stored as Python dictionary called pore_properties, thus the following will produce a list of all pore property lists currently stored in the dictionary:

>>> pn.pore_properties.keys()
['diameter', 'numbering', 'volume', 'seed', 'coords', 'Pc_invaded', 'type']

And similarly for throats, one gets the following:

>>> pn.throat_properties.keys()
['volume', 'diameter', 'numbering', 'connections', 'length', 'seed', 'Pc_invaded', 'Pc_entry', 'type']

A complete list of the default predefined pore and throat properties is given below:

**Common Pore and Throat Properties**

*'numbering'*: blah

*'type'*: blah

**Throat Specific Properties**

*'connections'*: blah

**Pore Specific Properties**

*'coords'*: blah

 





Examples
--------

To reserve space for a network with the default number of pores
and throats execute

>>> import OpenPNM as PNM
>>> net=PNM.Network.GenericNetwork()
>>> net.print_overview()
==================================================
= Overview of network properties
--------------------------------------------------
Basic properties of the network
- Number of pores:    10
- Number of throats:  20
Pore properties:
    numbering   int64     (10,)
Throat properties:
    connections int64     (20, 2)
    numbering   int64     (20,)
--------------------------------------------------

The following example plots the default adjacency matrix:

.. plot::
    
    import pylab as pl
    import OpenPNM
    net = OpenPNM.Generators.Cubic().generate()
    net.create_adjacency_matrix()
    pl.spy(net._adjmatrix)
    pl.show()



