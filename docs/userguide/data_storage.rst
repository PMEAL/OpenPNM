*******************************************************************************
Network Architecture and Data Storage Formats
*******************************************************************************

OpenPNM utilizes the object oriented capacities of Python.  The code is built upon the idea that a network is an object.  A network object contains both the data that describes the network properties, and the tools, functions, or methods needed to access this data in ways applicable to the pore network modeling paradigm.  One key feature of this object is that it is completely agnostic about the type of network it describes; a random, cubic or another netowrk topology is stored in exactly the same manner.  The most important repercussion of this choice is the fact that all physical algorithms (such as diffusion or drainage) operating on the network can be fully generic, and the all methods that read and write network property data can be fully generic as well.  

As the name suggests, pore network modeling borrows signifcantly from the fields of network and graph theory.  During the development of OpenPNM, it was debated whether existing Python graph theory packages (such as `graph-tool <http://graph-tool.skewed.de/>`_ and `NetworkX <http://networkx.github.io/>`_) should be used to store the network topology.  It was decided that storage of network property data could be handled very efficiently using 1D arrays, which allowed for a high degree of code vectorization.  Fortuitously, around the same time as this dicussion, Scipy started to include the 'compressed sparse graph' library, which contained numerous graph theory algorithms.  The CSGraph libary requires adjacency matrices in a compressed sparse storage scheme, which happens to be how OpenPNM stores network connections as descrbed below.

===============================================================================
Network Architecture
===============================================================================

-------------------------------------------------------------------------------
Pore and Throat Numbering
-------------------------------------------------------------------------------

One of the main design principle of OpenPNM was to accomodate pore networks of arbitrary dimensionality and shape.  Cubic networks are commonly used in pore network modeling, with each pore connected to 6 or 26 neighbors.  The type of network *can* be represented as cubic matrices in numerical simulations, and this has the advantage that it is easily interpreted by human users.  Representing networks this way, however, clearly lacks generality.  Networks extracted from tomographic images, or generated using random pore placements connected by Delaunnay tesselations require a different approach.  OpenPNM uses network representation schemes borrowed from graph theory, such as adjacency and incidence matrices, that can be used to represent *all* network topologies. 

The basic definitions used by OpenPNM are:

1. Pores can have an arbitrary number of throats, including zero

2. A throat connects exactly two pores, no more and no less

3. Two pores are connected by no more than one throat

A network has a certain number of pores, *Np*, and a certain number of throats, *Nt*.  Typically, *Nt* > *Np* since most pores have more than 1 throat.  If every pore has 1 throat (i.e. the network forms a linear chain), then *Nt* = *Np* - 1.  Of course, Nt can be zero but this would not be a useful network.  It can be *unofficially* stated that a network should have at least 2 pores connected by at least 1 throat (*Np* > 1 and *Nt* > 0).  

Each pore and throat in the network has a unique ID number.  In OpenPNM the ID number is *implied* by array element number, meaning that any information stored in element *i* of an array implicity applies to pore (or throat) *i*.  Or in other words, finding information about pore (or throat) 0 is accomplished by looking into element 0 of an array.  There is no correspondence between pore number and throat number, meaning that throat *i* may or may not connected with pore *i*.  Python uses 0-based array indexing so the ID numbers start at 0, which can be a source of confusion when representing connections using sparse representations,

-------------------------------------------------------------------------------
Adjacency Matrices
-------------------------------------------------------------------------------

When each pore has a unique ID number it is logical to store the network connectivity as a list pore ID numbers to which a given pore is connected.  Although it is possible to store these *lists* as a list-of-*lists*, graph theoreticians have devised an more elegant and powerful approach, which OpenPNM has adopted, called adjacency matrices.  An adjacency matrix is a 2D matrix of size *Np*-by-*Np*.  A value of 1 is placed at location (*i*, *j*) to indicate that pores *i* and *j* are connected.  In pore networks there is generally no difference between traversing from pore *i* to pore *j* or from pore *j* to pore *i*, so a 1 is also placed at location (*j*, *i*).  This means that determining which pores are connected directly to a given pore (say *i*) can be accomplished by finding the locations of non-zeros in row *i*.  In graph theory terminology this is deemed an *undirected* network, meaning that the *direction* of traversal is immaterial.  The adjacency matrix of an undirected network is symetric.  

Since the adjacency matrix is symetric it is redundant to store the entire matrix when only the upper triangular part is necessary.  Furthermore, because pores are generally only connected to nearby pores, the number of throats per pore is a very small faction of the total number of throats.  This means that there are very few non-zero elements on each row, so the adjacency matrix is highly sparse.  This fact naturally lends itself to sparse storage schemes.  OpenPNM uses uses the IJV sparse storage scheme to store the upper triangular poriton of the adjacency matrix.  The *IJV* scheme is simply an *Np*-by-3 array of the (*I*, *J*) coordinates of each nonzero element in the adjacency matrix, along with the corresponding non-zero value (*V*).  (The scipy.sparse module calls this the Coordinate or COO storage scheme, but it is more widely known as IJV).  For example, to denote a value of 1 on row 3 and column 7, the *IJV* storage scheme would include an entry IJV = [3, 7, 1].  Each nonzero element in the adjacency matrix corresponds to a row to the *IJV* array.  Moreover, the number of non-zeros in the upper triangular portion of the adjacency matrix is equal to the number of throats in the network, so the dimensions of the *IJV* array is *Nt*-by-3.  This is not a coincidence; a key feature of the adjacency matrix is that each non-zero element direclty corresponds to a throat.  Because throat numbers are implicitly defined by their location in an array, then the IJV sparse storage scheme automatically assigns throat ID numbers when the IJV array is generated.  For instance, when scanning the adjacency matrix from left-to-right, top-to-bottom, the first non-zero element encountered (say at location [0,5]) would be assigned throat number 0, and stored as IJV[0] = [0,5,1].  

One further optimization used by OpenPNM is to drop the V from the IJV format since the non-zeros in the adjacency matrix are all 1.  This results in a *Nt*-by-2 array which is called *connections*.  Any desired throat property array can be appended as a third column to the *connections* array to fully specify the IJV format for use with the scipy.sparse or scipy.csgraph functions.  OpenPNM provides a routine for this operation ('fill_adjacency_matrix'), which takes the desired throat property list to insert into *V* as an argument.  

In summary, when storing network connectivity as the upper triangular portion of an adjacency in the IJV sparse storage format, the end result is an *Nt*-by-2 list describing which pores are connected by a given throat.  These connections are a fundamental property associated with each throat in the same way as throat diameter or capillary entry pressure.  This highly distilled storage format minimized memory usage, allows for vectorization of the code, is the most efficient means of generating a sparse matrix, and corresponds perfectly with the storage of other throat properties using the ID number implicitly defined by the list element location. 

-------------------------------------------------------------------------------
Other Sparse Storage Schemes
-------------------------------------------------------------------------------
The IJV storage format corresponds perfectly with the way other throat data is stored in OpenPNM, however some tasks and queries are performed more efficiently using other storage formats.  OpenPNM converts between these formats internally as needed.  

-------------------------------------------------------------------------------
Incidence Matrices
-------------------------------------------------------------------------------
Another way to represent network connections is an incidence matrix.  This is similar to an adjacency matrix but rather than denoting which pores are connected to which, it denotes which pores are connected to which throats.  An incidence matrix is *Np*-by-*Nt* is size, with *Nt* nonzero elements.  The incidence matrix is useful for quickly querying which throats are connected to a given pore by finding the location of non-zero elements on a row.  Incidence matrices are generated as needed by OpenPNM internally for performing such queries, and the user does not usually interact with them.  

===============================================================================
Network Data Storage
===============================================================================
OpenPNM stores all pore and throat properties as Numpy ndarrays.  ndarrays are a numerical data type provided by the Numpy package (which is embedded in the Scipy package) that allow for the type of numerical manipulations that scientists and engineers expect, such as vectorization, slicing, boolean indexing and so on.  Pore properties are stored as arrays of size *Np*-by-*n), where *Np* is the number of pores in the network and *n* is almost always 1, (e.g. pore volume is stored as an *Np*-by-1 array), with a few expections (e.g. spatial coordinates are stored as *Np*-by-3 for 3-dimensional space).  Throat properties are almost always stored as *Nt*-by-*m* arrays where *Nt* is the number of throats in the network.  Again, *m* is almost always 1 with a notable exception being the connections property that is discussed in detail above. 

OpenPNM uses implied pore and throat numbering, meanin that the property for pore (or throat) *i* is stored in element *i* of the corresponding property array.  (Aside: It is conceivable that an extra array called 'index' could be used to remove the implicitness of the numbering.  For instance *pore_index* = [0,2,1,3] would indicate that the properties for 'pore 2' are located in element 1.  Given a list of pore properties pore_props = [0.1, 0.2, 0.3 0.4], one could extract an ordered list as *pore_props*[*index*] = [0.1, 0.3, 0.2, 0.4].  This extra layer of indexing is confusing and makes it more difficult to setup vectorized and boolean masked statements.)

To examine the properties of a network, start by generating a small network of 3-by-3-by-3 as follows:

.. code-block:: python
   import OpenPNM
   pn = OpenPNM.Geometry.Cubic(divisions=[3,3,3]).generate()

This creates a cubic network with 27 pores and 54 throats.  A quick summary of the network data can be had by typing:

>>> print pn

The following output will be produced:

.. code-block:: python

    ==================================================
    Overview of network properties
    --------------------------------------------------
    Basic properties of the network
    - Number of pores:   27
    - Number of throats: 54

    Pore properties:
        diameter            float64             (27L,)              
        numbering           int32               (27L,)              
        volume              float64             (27L,)              
        seed                float64             (27L,)              
        coords              float64             (27L, 3L)           
        type                int8                (27L,)              
    Throat properties:
        volume              float64             (54L,)              
        diameter            float64             (54L,)              
        numbering           int64               (54L,)              
        connections         int32               (54L, 2L)           
        length              float64             (54L,)              
        seed                float64             (54L,)              
        type                int8                (54L,) 
        
As can be seen, the network generation produces several basic pore and throat properties by default.  Note that the length of the pore and throat property lists correspond to the number of pores and throats in the network (27 and 54 respectively).  Most of the data are stored in 1D arrays, with two exceptions.  The pore property 'coords' gives the spatial location of the pore center in 3D cartesian coordinates, so each pore requires a set of X, Y and Z values.  The throat property 'connections' gives the ID number of the two pores it connects, or in other words it gives the IJ portion of the IJV sparse storage of the adjacency matrix.  

These data arrays are stored as part of the network object using Python dictionaries.  A Python dictionary is a form of structured variable where each entry in the dictionary has a 'key' : {value} pair.  The 'key' is the name of the of the value, and the {value} can be any data type.  In OpenPNM the values are all ndarrays.  

The dictionaries are called pore_properties and throat_properties.  To access the diameter of pores type:

>>> pn.pore_properties['diameter']

And similarly for throats:

>>> pn.throat_properties['diameter']

A quick way to find all properties currently stored in a dictionary is the .keys() method as follows:

>>> pn.pore_properties.keys()
['diameter', 'numbering', 'volume', 'seed', 'coords', 'type']

-------------------------------------------------------------------------------
Mandatory Pore & Throat Properties
-------------------------------------------------------------------------------
The default setup of the Cubic generator produces a number of pore and throat properties based on commonly used assumptions.  Only a few of these properties are truly essential to defining the pore network.  

**'connections' and 'coords'**

The spatial position of each pore is obviously a defining feature of a given pore network, so the 'coords' pore property is essential.  Equally essential to defining a network is the 'connections' throat property since this describes how the pores are connected or networked.  From a physical point of view, these are the only properties required to define a basic (though not very functional network).  With this informaiton it would be possible to generate a 3D images of the pore and throat network.  

**'type' and 'numbering'**

The 'type' and 'numbering' properties are also considered mandatory since OpenPNM relies on these for various internal calculations and network queries.  

The 'numbering' array is actually somewhat redundant since pore and throat numbers are implicitly defined by their array location.  This array is quite useful for boolean mask logic to find pores that meet a specific criteria.  For instance, to find all pores whose diameter is below average type:

>>> dia_mean = sp.mean(pn.pore_properties['diameter'])
>>> mask = pn.pore_properties['diameter'] < dia_mean
>>> small_pores = pn.pore_properties['numbering'][mask]
>>> small_pores
array([ 0,  3,  5,  6,  7,  9, 10, 12, 13, 14, 16, 17, 20, 21, 22])

(Note that the pore diameters are assigned randomly, so different network realizations will have different 'small_pores')

The 'type' property is used by OpenPNM to differentiate between internal pores and boundary pores (and throats).  A 'type' value of zero indicates an internal pore, and a value > 0 indicates a boundary pore.  Boundary pores are further distinguished by values between 1 and 6 to indicate on which boundary they lie: 1 and 6 for z-faces, 2 & 5 for x faces and 3 & 4 for y faces.  This convention was inspired by the number on dice, where opposite sides all add up to 7.  Obviously, this numbering boundary pores in this way implies a cubic network domain, which may not always be the case.  

-------------------------------------------------------------------------------
Common Pore & Throat Properties
-------------------------------------------------------------------------------
The GenericGeometry class includes several methods that produce some additional pore and throat properties beyond the mandatory ones described above.  These including this like 'diameter' and 'volume'.  


-------------------------------------------------------------------------------
Adding New Pore & Throat Properties
-------------------------------------------------------------------------------



===============================================================================
Querying Network Data and Properties
===============================================================================

The OpenPNM network object not only stores the network data, but also contains numerous methods for extracting information about the network from that data.  


.. automethod:: OpenPNM.Network.GenericNetwork.get_num_pores()

.. automethod:: OpenPNM.Network.GenericNetwork.get_num_throats()

.. automethod:: OpenPNM.Network.GenericNetwork.get_neighbor_pores()

.. automethod:: OpenPNM.Network.GenericNetwork.get_neighbor_throats()

.. automethod:: OpenPNM.Network.GenericNetwork.get_num_neighbors()

.. automethod:: OpenPNM.Network.GenericNetwork.get_connected_pores()

.. automethod:: OpenPNM.Network.GenericNetwork.get_connecting_throat()




