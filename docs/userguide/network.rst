.. _network:

###############################################################################
Network
###############################################################################
This module is the heart of OpenPNM.  It contains the ``GenericNetwork`` class which possesses a suite of network query methods, based the graph theory concepts of adjacency and incidence matrices.  The methods in ``GenericNetwork`` are fully agnostic to the type and topology of network due the generalized way that OpenPNM stores data.  This is explained in more detail in the :ref:`here<data_storage>`.

The ``GenericNetwork`` class on its own has no topology.  If you instantiate a ``GenericNetwork`` it will have no pores or throats:

>>> pn = OpenPNM.Network.GenericNetwork()
>>> pn.num_pores()
0
>>> pn.num_throats()
0

You can get a quick overview of the network properties by 'printing' it on the command line:

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

As can be seen, a basic empty network has 0 pore coordinates and 0 throat connections, and the label 'all' exists but is applied nowhere.  

The network module contains numerous subclasses of ``GenericNetwork``, which possess the code for actually generating specific network topologies (e.g. cubic, random, etc).  All subclasses derive from ``GenericNetwork`` so have its methods, as well as any additional methods relevant to the specific topology.  Generating a standard Cubic network is accomplished with:

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

-------------------------------------------------------------------------------
Properties
-------------------------------------------------------------------------------
Accessing the pore and throat property can be accomplished using standard Python ``dict`` syntax:

>>> pn['pore.index']
array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
       17, 18, 19, 20, 21, 22, 23, 24, 25, 26])
>>> pn['pore.index'][0]
0
>>> pn['pore.index'][[0,1,2]]
array([0, 1, 2])

It's also possible to write to the object:

>>> pn['pore.index'] = 1
>>> pn['pore.index']
array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1])
>>> pn['pore.blah'] = 2

This statement reveals an important aspect of the OpenPNM framework that was introduced above.  All 'pore' arrays are forced to be Np long, so in the absence of specifically indexing into the array the entire array is set to the scalar value.  Specific pores can be set with:

>>> pn['pore.blah'][[0,1,2,3]] = 0
>>> pn['pore.blah']
array([0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
       2, 2, 2, 2])

All the main OpenPNM objects have a method for quickly listing all of the defined pore and throat properties using ``props``.  With no arguments this returns all properties, or it can return just pore or throat properties:

>>> pn.props()
['pore.blah', 'throat.conns', 'pore.index', 'pore.coords']
>>> pn.props('pore')
['pore.blah', 'pore.index', 'pore.coords']
>>> pn.props('throat')
['throat.conns']

This is useful for iterating through all properties on the object, or just for visually inspecting the object.

-------------------------------------------------------------------------------
Labels
-------------------------------------------------------------------------------
The print-out of ``Cubic`` also includes a number of labels that were automatically applied by the generator. Labels are quite useful as they allow a quick way to select a subset of pores:

>>> pn.pores('pore.back')
array([18, 19, 20, 21, 22, 23, 24, 25, 26], dtype=int64)
>>> pn.pores(['pore.back','pore.front'])
array([ 0,  1,  2,  3,  4,  5,  6,  7,  8, 18, 19, 20, 21, 22, 23, 24, 25,26], dtype=int64)

Note that this could also have been achieved by checking pore coordinates and filtering based on their location, which is how the generator applies the labels initially.  Any complicated query used to find pores or throats can be stored as a label for future use:

>>> Ps = pn['pore.coords'][:,2] > sp.mean(pn['pore.coords'][:,2])
>>> pn['pore.top_half'] = Ps
>>> pn.pores('pore.top_half')
array([ 2,  5,  8, 11, 14, 17, 20, 23, 26], dtype=int64)

-------------------------------------------------------------------------------
Topology Queries
-------------------------------------------------------------------------------
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

-------------------------------------------------------------------------------
Topology Manipulations and Operations
-------------------------------------------------------------------------------
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