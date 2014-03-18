.. _network:

###############################################################################
The Network Object
###############################################################################

The **Network** is the central object in OpenPNM.  This object contains 3 main aspects:

1.  **Query methods**:  There are methods for simple queries like finding the number of pores or throats in the network.  There are also more complex methods for finding the number of neighbor pores connected to a given pore or set of pores.  

2.  **Adjacency and Incidence Matrices**:  These matrices are stored on the network object as Scipy sparse arrays. They can be created using the ``create_adjacency`` or ``create_indicidence`` methods.  They can be accessed directly for use in algorithms as needed (as in ordinary percolation).  Some of the *query* type methods access these matrices to find connected throats and neighbors.  

3.  **Data Storage**:  All of the data pertaining to topology and geometry are stored on the Network object under ``_pore_data`` and ``_throat_data``.  These dictionaries are filled with Numpy ndarrays, each stored under a descriptive *key* such as 'diameter' or 'volume'.  These arrays should not be accessed directly however.  The **Network** inherits the *setter* and *getter* methods from **Tools**, so these should be used.