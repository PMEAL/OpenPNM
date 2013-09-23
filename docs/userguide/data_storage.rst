===============================================================================
Network Architecture and Data Storage Formats
===============================================================================
OpenPNM utilizes the object oriented capacities of Python.  The code is built upon the idea that a network is an object.  A network object contains both the data that describes the network properties, and the tools, functions, or methods needed to access this data in ways applicable to the pore network modeling paradigm.  One key feature of this object is that it is completely agnostic about the type of network it describes; a random, cubic or another netowrk topology is stored in exactly the same manner.  The most important repercussion of this choice is the fact that all physical algorithms (such as diffusion or drainage) operating on the netowrk can be fully generic, and the fact that all methods that read and write network property data can be fully generic as well.  

As the name suggests, pore network modeling borrows signifcantly from the fields of network and graph theory.  During the development of OpenPNM, it was considered whether existing network tools (such as graph-tools and networkX) should be used to store the network topology.  It was decided that storage of network property data could be handled very efficiently using 1D arrays, which allowed for a high degree of code vectorization.  Fortuitously, Scipy released a version that contained the 'compressed sparse graph' library, which contained numerous graph algorithms.  The CSGraph libary requires adjacency matrices in a compressed sparse storage scheme, which happens to be how OpenPNM stores network connections.  

-------------------------------------------------------------------------------
Network Data Storage
-------------------------------------------------------------------------------
OpenPNM stores all pore and throat properties as Numpy ndarrays (a Numpy numerical data type).  Pore properties are stored as arrays of size Np x n, where Np is the number of pores in the network and n is almost always 1, (e.g. pore volume is stored as an Np x 1 array), with a few expections (e.g. spatial coordinates are stored as Np x 3 for 3-dimensional space).  Throat properties are almost always stored as Nt x m arrays where Nt is the number of throats in the network.  Again, m is almost always 1 with a notable exception being the connections property that is discussed in detail below. 

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

*******************************************************************************
Network Storage
*******************************************************************************






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



