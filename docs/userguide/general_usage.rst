.. _general:

===============================================================================
General Overview
===============================================================================
Before diving into the detailed workings of each Core object, it is worthwhile to give an overview of the general behavior of these objects.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The Core Class
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The Core objects in an OpenPNM simulation descend from the Python ``dict`` class, and have additional methods added.  All of the main OpenPNM objects are children of Core, and hence are often referred to collectively as 'Core' objects.  This section will outline features common to the Core objects.  

-------------------------------------------------------------------------------
Data Storage
-------------------------------------------------------------------------------
Each OpenPNM Core object is a ``dictionary`` which is the Python equivalent to a structured variable or struct in other languages.  This allows data to be stored on the object and accessed by name, with a syntax like ``network[‘pore.diameter’]``.  All pore and throat data are stored as 1D vectors or lists, of either Np or Nt length representing the number of pores and throats in the Network (or object), respectively.  This means that each pore (or throat) has a number, and all properties for that pore (or throat) are stored in the array corresponding to that number.  Thus, the diameter for pore 15 is stored in the ‘pore.diameter’ array in element 15, and the length of throat 32 is stored in the ‘throat.length’ array at element 32.  All data is converted to a Numpy array when stored in the dictionary.  Numpy is the *de facto* standard numerical data type in Python.  Not only does this list-based approach allow for any topology to be stored in an equivalent manner, but it is optimal for vectorised coding, thus enabling high performance numerical operations with the Numpy and Scipy libraries which support vectorised, element-wise operations on multicore processors.   

Several rules have been implemented to control the integrity of the data, as described in the note below.  Firstly, all list names must begin with either ‘pore’ or ‘throat’ which obviously serves to identify the type of information stored there.  Secondly, for the sake of consistency only data arrays of length Np or Nt are allowed in the dictionary.  This rule forces the user to be cognizant of the list-based numbering scheme used to identify pores and throats.  Any scalar values written to the dictionaries are broadcast into full length vectors, effectively applying the scalar value to all locations.  Attempts to write only partial data to a subset of pores or throats will result in the assignment of ‘NaN’ values (Not-a-Number) to other locations.  Finally, any data that is Boolean will be treated as a ``label`` while all other numerical data is treated as a ``property``.  The difference between these is outlined below.  

Data is stored in a compartmentalized fashion, with topological information on the Network object, and geometrical information the Geometry objects.  Though very helpful, this practice leads to a possible point of confusion: since there can be and usually are multiple Geometry objects the number of pores (and throats) stored on each Geometry object differs from the total number of pores (and throats) in the full Network.  Moreover, on each individual Core object the pores are stored using their own internal numbering, so pore 100 on the Network may be pore 1 on a Geometry object.  There are methods for mapping between objects (``map_pores`` and ``map_throats``) to help track these numbering mismatches.  

.. note:: **__setitem**

    ``__setitem_`` is the private method on ``dict`` that is called when the dictionary syntax is used to write values, so ``pn['pore.test'] = 0`` is equivalent to ``pn.__setitem__('pore.test',0)``.  OpenPNM subclasses the ``__setitem__`` method to intercept data and ensure it meets certain criteria before being written to the objects.  The two main rules are that (1) all dictionary keys must start with either 'pore' or 'throat', and (2) all data must be of the correct length, either Np or Nt long, where Nt is the number of throats and Np is the number of pores on the object.

-------------------------------------------------------------------------------
Properties and Labels
-------------------------------------------------------------------------------
OpenPNM differentiates between two types of data for pores and throats: 'properties' and 'labels'.  The only difference between these is that 'labels' are Boolean arrays (True / False), while 'properties' are numerical data types.  

The physical details about pores and throats is referred to as 'properties', which includes information such as 'pore volume' and 'throat length'.  Properties can be accessed using standard Python dictionary syntax:

>>> pn['pore.index'][1]
1
>>> pn['pore.index'][[0,1]]
array([0, 1])

Writing data also uses dictionary syntax, but with a few caveats due to the fact that OpenPNM has subclassed ``__setitem__`` to 'protect' the integrity of the data. 

>>> pn['pore.index'][10] = 3
>>> pn['pore.index'][10]
3

The main 'caveat' is that data will all be forced to be either Np or Nt long, so the following attempt to write a scalar value will result in a vector of length Np (filled with 1's): 

>>> pn['pore.dummy'] = 1.0

To quickly see a list of all defined 'properties' use ``props``.  You can specify whether only 'pore' or 'throat' properties should be returned, but the default is both:

>>> pn.props()
['pore.index', 'pore.coords', 'throat.conns']
>>> pn.props('pore')
['pore.index', 'pore.coords']

For more details on ``props``, see the method's docstring.  

The second type of information is referred to as 'labels'.  Labels were conceived as a means to dynamically create groups of pores and throats so they could be quickly accessed by the user.  For instance, in a Cubic Network it is helpful to know which pores are on the 'top' surface.  This label is automatically added by the topology generator, so a list of all pores on the 'top' can be retrieved by simply querying which pores possess the label 'top'.  

The only distinction between 'labels' and 'properties' is that 'labels' are boolean masks of True/False.  Thus a True in element 10 of the array 'pore.top' means that the label 'top' has been applied to pore 10.  Adding and removing existing labels to pores and throats is simply a matter of setting the element to True or False.  Creating a new label is a bit more tricky.  'label' arrays are like any array and they must be defined before they can be indexed, so to apply the label 'dummy_1' to pore 10 requires the following 2 steps:

>>> pn['pore.dummy_1'] = False
>>> pn['pore.dummy_1'][10] = True

Now that this label array has been created and True values have been inserted, it is a simple matter to recall which pores have 'dummy_1' by finding the locations of the True elements:

>>> sp.where(pn['pore.dummy_1'])[0]

OpenPNM provides a more convenient way to perform this query with the ``pores`` and ``throats`` methods that are outlined below.  

The ``labels`` method can be used to obtain a list of all defined labels. This method optionally accepts a list of pores or throats as an argument and returns only the labels that have been applied to the specified locations.  

>>> pn.labels()
['pore.all', 'pore.back', 'pore.bottom', 'pore.front', 'pore.internal', 'pore.left', 'pore.right', 'pore.top', 'throat.all']

``labels`` also has a ``mode`` argument that controls some set-theory logic to the returned list (such as 'union', 'intersection', etc).  See the method's docstring for full details.

-------------------------------------------------------------------------------
Counts and Indices
-------------------------------------------------------------------------------
One of the most common questions about a network is "how many pores and throats does it have?"  This can be answered very easily with the ``num_pores`` and ``num_throats`` methods.  Because these methods are used so often, there are also shortcuts: ``Np`` and ``Nt``.  

>>> pn.num_pores()
27

It is also possible to 'count' only pores that have a certain label (shortcuts``Np`` and ``Nt`` don't work with this counting method):

>>> pn.num_pores('top')
9

These counting methods actually work by counting the number of True elements in the given label array.  

Another highly used feature is to retrieve a list of pores or throats that have a certain label applied to them.  This is of course the entire purpose of 'labels'.  To receive a list of pores on the 'top' of the Cubic network:

>>> pn.pores('top')
array([ 2,  5,  8, 11, 14, 17, 20, 23, 26], dtype=int64)

The ``pores`` and ``throats`` methods both accept a 'mode' argument that allows for set-theory logic to be applied to the query, such as returning 'unions' and 'intersections' of locations. For complete details see the docstring for these methods.  

Often, one wants a list of all pore or throat indices on an object, so there are shortcut methods for this: ``Ps`` and ``Ts``.

.. note:: **The Importance of the 'all' Label**

   All objects are instantiated with a 'pore.all' and a 'throat.all' label.  These arrays are essential to the framework since they are used to define how long the 'pore' and 'throat' data arrays must be.  In other words, the ``__setitem__`` method checks to make sure that any 'pore' array it receives has the same length as 'pore.all'.  Moreover, the ``pores``, ``throats``, ``num_pores`` and ``num_throats`` methods all have the label 'all' as their default so if no label is sent 'all' pores or throats are considered.  

-------------------------------------------------------------------------------
Naming
-------------------------------------------------------------------------------
All OpenPNM objects are given a name upon instantiation.  The name can be specified in the initialization statement:

>>> pn = OpenPNM.Network.Cubic(shape=[3,3,3],name='test_net_1')
>>> pn.name
'test_net_1'

The name of an object is stored under the attribute 'name'. If a name is not provided, then a name will be automatically generated by appending 5 random characters to the class name (e.g. 'Cubic_riTSw').  It is also not possible to have two objects with the same name associated with a Network.  Names can be changed by simply assigning a new string to the ``name``.  
   
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Models Dictionary
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Models are one of the most important aspects of OpenPNM, as they allow the user to specify a 'model' for calculating 'pore.volume', rather than just entering numerical values into a geometry_object['pore.volume'] array.  It is mainly through customized models that users can tailor OpenPNM to a specific situation, though OpenPNM includes a variety of pre-written models.  These are stored under each Module in a folder called 'models'.  For instance, Geometry.models.pore_diameter contains several methods for calculating pore diameters.  

Each Core object has a ``models`` attribute where all information about pore-scale models are stored.  Upon instantiation of each ``Core`` object, a ``ModelsDict`` object is stored in its ``models`` attribute.  The ``ModelsDict`` class is a subclass of the Python ``dict`` class, which has several features added for dealing specifically with models.  A detailed description of the Models Dictionary class can be found :ref:`here<models>`.

Adding a model to an object 

The ``add_model`` method accepts 3 main types of argument.  

(1) It needs to know which pore or throat property does this model calculate.  This is the *propname* argument, and would be 'pore.seed' or 'pore.diameter' for the example above.
(2) It needs the actual function that should be used.  In Python it is possible to pass a function as an argument as easily as passing an integer.  The  *model* argument should be a handle to the function of choice such as Geometry.models.pore_size.sphere.
(3) It can optionally accept an arbitrary number of arguments that will be passed directly to the 'model'.  

These 3 requirements are well demonstrated by the random pore seed model:

.. code-block:: python

	geom = OpenPNM.Geometry.GenericGeometry()  # Creates an empty Geometry object
	mod = OpenPNM.Geometry.models.pore_misc.random  # Get a handle to the desired model
	geom.add_model(propname='pore.seed',model=mod,seed=0)  # Assign model to the object
	
The *propname* and *model* arguments are required by the ``add_model`` method, but the *seed* argument is passed on the model, and it specifies the initialization value for the random number generator.  

The ``add_model`` method actually runs the model and places the data in the dictionary given by *propname*. It also saves the model in a special dictionary attached tyo the object (object.models) also under the same *propname*.  When the data is requested from the object it returns the 'static' copy located in the object's dictionary.  In order to recalculate the data the model stored in the private dictionary must be rerun.  This is accomplished with the ``regenerate`` method.  This method takes an optional list of *propnames* that should be regenerated.  It should also be pointed out that models are regenerated in the order that they were added to the object so some care must be taken to ensure that changes in property values cascade through the object correctly.  



















