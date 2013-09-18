===============================================================================
Introduction
===============================================================================
OpenPNM is a pore network modeling package coded in the Python programming language.  OpenPNM makes extensive use of the Scipy and Numpy libraries, which add considerable scientific and numerical computation capabilities to Python.  

-------------------------------------------------------------------------------
OpenPNM Architecture
-------------------------------------------------------------------------------
OpenPNM utilizes the object oriented capacities of Python.  The code is built upon the idea that a network is an object.  A network object contains both the data that describes the network properties, and the tools, functions, or methods needed to access this data in ways applicable to the pore network modeling paradigm.  One key feature of this object is that it is completely agnostic about the type of network it describes; a random, cubic or another netowrk topology is stored in exactly the same manner.  The most important repercussion of this choice is the fact that all physical algorithms (such as diffusion or drainage) operating on the netowrk can be fully generic, and the fact that all methods that read and write network property data can be fully generic as well.  

As the name suggests, pore network modeling borrows signifcantly from the fiels of network and graph theory.  During the development of OpenPNM, it was considered whether existing network tools (such as graph-tools) should be used to store the network topology.  It was decided that storage of network property data could be handled very efficiently using 1D arrays, which allowed for a high degree of code vectorization.  Fortuitously, Scipy released a version that contained the 'compressed sparse graph' library, which contained numerous graph algorithms.  The CSGraph libary requires adjacency matrices in a compressed sparse storage scheme, which happens to be how OpenPNM stores network connections.  

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

Most of these are self explanatory, but a few are deserving of careful explanation.  

**'numbering'**
The *'numbering'* property in both pore and throat arrays is simply an enumerated list from 0 to Np or Nt.  These numbers correspond to the pore or throat number, which may seem redundant this also happens to to corresond the location in the array.  Numbering is useful in some cases where boolean arrays are used to select a subset of pores (or throats).  The boolean mask can be input into the numbering property to produce a list of pore (or throat) numbers where the condition is true.  

**'coords'**
The spatial location of each pore body is given by this *'coords'* property.  Each element in this array is 3 columns wide, containing the x, y and z components of the pore location in Cartesean coordinates.  

**'seed'**
This seed is used in the specified statistical distribution to obtain pore (and throat) sizes.  In the simplest case, seed is selected randomly for each pore, and each throat receives the smaller of it's neighbors seeds.  This property can be generated in more complex manner to include spatial correlations or other features.  

**'type'**
The *'type'* property is 
