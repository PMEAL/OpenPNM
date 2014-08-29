.. _getting_started:

###############################################################################
Getting Started with OpenPNM
###############################################################################
The OpenPNM framework is build upon 5 main objects, each outlined below.  All of these objects are a subclass of the Python 'dictionary' (``dict``), which is a data storage class similar to 'structs' in C or Matlab.  This essentially means that multiple pieces of data can be stored on each object, and accessed by name (i.e. object_name['pore.diameter']).  Using ``dict`` as a starting point, the OpenPNM classes provide an easy to use and understand data storage scheme that allows direct access and interaction with the network data.  The OpenPNM framework augments the basic ``dict`` by adding a variety of methods that work with the data stored in the dictionaries.  For instance, the Python ``dict`` class only has a few basic methods for reading and removing data from the dictionary, while an OpenPNM object has many additional methods related to pore network modeling, such as querying the number of pores in the network.  

Before proceeding with an overview of the framework, a few words must be said about how OpenPNM stores data, the type of data it stores, and some basic rules it applies to data formats and names.  Firstly, all pore and throat data are stored as lists (actually Numpy arrays).  If a network has 100 pores then all pore data (i.e. diameter) is stored in a list 100 elements long.  Each pore and throat is associated with a number which corresponds to the location in the pore or throat lists where its information is stored, so data for pore 10 is stored in element 10.  This list format allows any topology to be described by the framework, and also enables easy vectorization of code for performance gains.  A detailed description of this storage scheme is given :ref:`here <data_storage>`.  

Secondly, OpenPNM differentiates between two types of data for pores and throats.  The physical details about pores and throats is referred as *properties*, and this includes information such pore volume and throat length.  The second type of information is referred to as *labels*.  Labels were conceived as a means to dynamically create groups of pores and throats so they could be quickly accessed by the user.  For instance, in a cubic network it is helpful to know which pores are on the 'top' surface.  This label is automatically added by the topology generator, so a list of all pores on the 'top' can be retrieved by simply querying which pores possess the label 'top'.  

Finally, OpenPNM forces all data and label *names* to be either of the form 'pore.name' or 'throat.name'.  All 'pore.name' lists are Np long, and all 'throat.name' lists are Nt long, where Np and Nt are the number of pores and throats.  Forcing all arrays to be of the same length ensures that vectorized code operations will always work, and enforcing this naming allows the framework to make each array to correct length.  It is also very convenient to see what type of data is stored in a list just by glancing at its name.  

===============================================================================
Main Modules
===============================================================================

1 `Network`_: The Network object is the main *controller* object for the framework.  It contains all the topological information about the **Network**, as well as methods for querying and manipulating the topology. 

3 `Geometry`_: Geometry objects control the pore-scale geometrical properties of the network such as pore volume and throat diameter.  A simulation may have multiple Geometry objects depending on the problem being modeled.  For instance, a stratified material may have a separate Geometry object for each layer if the pore and throat sizes differ between them.  

4 `Phases`_: Phase objects contain information about the thermo-physical properties of the liquids, gases, solids that exist in the pore space.  For instance, a Phase object for water would possess its temperature, as well as models for calculating its viscosity as a function of temperature (and any other relevant properties).

5 `Physics`_: Physics objects contain methods for calculating pore scale physics properties which combine fluid and geometry values.  For instance, the hydraulic conductance of a throat requires knowing the throat diameter and length, as well as the fluid viscosity.  

6 `Algorithms`_: This module is the home of the actual algorithms that use the network properties defined by the above modules.  OpenPNM ships with a good assortment of standard algorithms, but is meant to be extended by users adding custom algorithms.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Network
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This module is the heart of OpenPNM.  It contains the ``GenericNetwork`` class which possesses a suite of network query methods, based the graph theory concepts of adjacency and incidence matrices.  The methods in ``GenericNetwork`` are fully agnostic to the type and topology of network due the generalized way that OpenPNM stores data.  This is explained in more detail in the :ref:`here<data_storage>`.

The ``GenericNetwork`` class on its own has no topology.  If you instantiate a ``GenericNetwork`` it will have no pores or throats:

>>> pn = OpenPNM.Network.GenericNetwork()
>>> pn.num_pores()
0
>>> pn.num_throats()
0

You can get a quick overview of the network properties by 'printing' it:

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


The network module contains numerous subclasses of ``GenericNetwork``, which possess the code for actually generating specific network topologies (e.g. cubic, random, etc).  All subclasses derive from ``GenericNetwork`` so have its methods, as well as any additional methods relevant to the specific topology.  Generating a simple cubic network is accomplished with:

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

The print-out of the network information shows that it has 27 pores and 54 throats, with properties of 'pore.coords', 'pore.index' and 'throat.conns'.  Because the ``Cubic`` class only generates the topology there is not any information about pores and throat sizes.  The only requirements of a topology are that the pores have spatial locations (given by 'pore.coords') and throats know which two pores they connect ('throat.conns').  ('pore.index' is required for other purposes).  

The print-out also includes a number of labels that were automatically applied by the generator. Labels are quite useful as they allow a quick way to select a subset of pores:

>>> pn.pores('pore.back')
array([18, 19, 20, 21, 22, 23, 24, 25, 26], dtype=int64)
>>> pn.pores(['pore.back','pore.front'])
array([ 0,  1,  2,  3,  4,  5,  6,  7,  8, 18, 19, 20, 21, 22, 23, 24, 25,26], dtype=int64)

Note that this could also have been achieved by checking pore coordinates and filtering based on their location, which is how the generator applies the labels initially.  Any complicated query used to find pores or throats can be stored as a label for future use:

>>> Ps = pn['pore.coords'][:,2] > sp.mean(pn['pore.coords'][:,2])
>>> pn['pore.top_half'] = Ps
>>> pn.pores('pore.top_half')
array([ 2,  5,  8, 11, 14, 17, 20, 23, 26], dtype=int64)

Network objects have many other methods for querying the topology, such as finding the neighbors of a pores, or finding the throat that connects 2 pores:

>>> pn.find_neighbor_pores(pores=[0])
array([1, 3, 9])
>>> pn.find_connecting_throat(P1=[0,0,0],P2=[1,3,9])
[[0], [18], [36]]

The best way to explore the available methods is to use an IDE editor that support the autocomplete function, such as Spyder.  This way, you can type ``pn.`` and a pop-up list of available methods will appear.  Extensive documentation is also included inside the OpenPNM code itself in the form of 'docstrings' which will be interpreted by Spyder and shown in the *Object Inspector*.  These docstrings give a description of the required and optional arguments to each method, along with examples and notes where applicable.  

More details on the Network object, such as how to make your own subclassed topology, are given :ref:`here<network>`.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Geometry
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The *Geometry* module controls all the network pore and throat size information.  This module contains the ``GenericGeometry`` class, which like all OpenPNM objects is subclass of Python's ``dict`` class, but has numerous OpenPNM specific methods added to it.  An empty ``GenericGeometry`` object *can* be initialized with no arguments, but this is not a useful object since it isn't associated with a network or assigned to any pores.  A more useful *Geometry* object would is obtained by instantiating a non-empty network, then assign a GenericGeometry to all pores and throats:

>>> pn = OpenPNM.Network.Cubic(shape=[3,3,3])
>>> Ps = pn.pores('pore.all')
>>> Ts = pn.throats('throat.all')
>>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,pores=Ps,throats=Ts)
>>> print(geom)
------------------------------------------------------------
OpenPNM.Geometry.GenericGeometry: 	GenericGeometry_ZpKsC
------------------------------------------------------------
#     Properties                          Valid Values
------------------------------------------------------------
1     pore.map                               27 / 27   
2     throat.map                             54 / 54   
------------------------------------------------------------
#     Labels                              Assigned Locations
------------------------------------------------------------
1     pore.all                            27        
2     throat.all                          54        
------------------------------------------------------------

Note that, because all objects in OpenPNM descend from the same base class, they all have many of the same query methods:

>>> geom.num_pores()
27
>>> geom.num_throats()
54
>>> geom.labels()
['pore.all', 'throat.all']
>>> geom.props()
['throat.map', 'pore.map']

This Geometry object is now associated with the network, ``pn``, and applied to all the pores and throats in the network.  At this point, however, it is still an empty object with no pore or throat size *information*.  To begin assigning size information, it is possible to simply assign values:

>>> geom['pore.seed'] = sp.rand(geom.Np)
>>> geom['throat.constant'] = 1.4

There are, however, very few cases where such a simple assignment is sufficient and in most cases, more elaborate pore scale models will be invoked.  For instance, it is common for throats to adopt the smaller of the seed values in it's two neighboring pores.  OpenPNM includes a library of pre-written models.  Pore scale geometry models are located under ``OpenPNM.Geometry.models``.  There are numerous files in this library with names that indicate their contents (i.e. pore_volume), and each of these files contain a variety of functions for calculating that property.  Specifying which models to use for a given property is done using the ``add_model`` method:

>>> import OpenPNM.Geometry.models as gm
>>> geom.add_model(propname='pore.seed',model=gm.pore_misc.random)

The above call to ``add_model`` does several important things.  Firstly, it creates a dictionary on ``geom`` called 'pore.seed'.  Secondly, it runs the function it received for the model argument and stores the returned values in 'pore.seed'.  Finally, it saves the model in a private dictionary on the `geom` object.  This final step is essential so that the *Geometry* object can retain a memory of it's models.  This means that the geometry can be *regenerated* using ``regenerate`` to recalculate its properties, and it means the object can be saved to disk and will still function fully when it's reloaded.

There can be multiple *Geometry* objects defined for different locations in a *Network* simultaneously.  This was intended to allow for multi-layer media (such fuel cell gas diffusion layers with micro-porous layers on one side), but is also quite useful when applying boundary pores which usually need to have special pore geometry such as 0 volume to produce consistent results.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Phases
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The **Phases** class contains methods for estimating or predicting the thermo-physical properties of fluids, such as viscosity or density, and this class can also be used to calculate solid properties.  (The name 'phase' was chosen over 'fluid' to be more general).

**Phase** objects are 'built' by the user to contain the specific methods that are to be used to calculate the phase properties.  For instance, a **Phase** object can calculate viscosity assuming a constant value, or Reynolds equation.  The user can use the methods supplied with OpenPNM, or add their own.  

Typically there will be multiple **PHase** objects defined for each simulation, since most models will have at least an invading fluid and a defending fluid.  There can be an unlimited number of phases associated with a **Network**.  

This module works identically to the `Geometry`_ module.    **Phase** objects are essentially empty when initialized, and the user adds the desired methods.  For example, the molar density of a fluid is a fundamental property that is required for many other calculations.  OpenPNM comes with several models that can calculate the molar density of a fluid.  These stored in a submodule called **molar_density**, and options include ``ideal_gas`` and ``real_gas``.  It is envisioned that each user may add their own methods for generating seeds (or ay other calculation) to this submodule (or any other), so the process of adding methods to the object was designed to be as flexible and customizable as possible.  For more information on this process see the :ref:`fluids <Fluids>` documentation.  The ``regenerate`` method does as its name suggests and regenerates all the data of the object.  As methods are added to the object (using ``add_method``) they are noted in a list, which is then referenced at regeneration time so all custom added methods are invoked. 

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Physics
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
asdf

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Algorithms
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
asdf

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Visualization
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
asdf

