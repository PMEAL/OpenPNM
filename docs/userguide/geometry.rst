.. _geometry:

===============================================================================
Geometry
===============================================================================
The *Geometry* module manages the network pore and throat size information.  This module contains the ``GenericGeometry`` class, which like all OpenPNM objects is subclass of Python's ``dict`` class, but has numerous OpenPNM specific methods added to it.  

.. inheritance-diagram:: OpenPNM.Geometry.GenericGeometry

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Basic Usage
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
An empty ``GenericGeometry`` object *can* be initialized with no arguments, but this is not a useful object since it isn't associated with a network or assigned to any pores.  A more useful *Geometry* object is obtained by instantiating a non-empty network, then assign a GenericGeometry to *all* pores and throats:

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

There are a few important steps that occur upon instantiation, aside from those outlined in the :ref:`General Overview<general>`.  Firstly, the Geometry stores a list of the Network pores and throats to which it applies.  These are stored under 'pore.map' and 'throat.map'.  For instance:

>>> pn = OpenPNM.Network.Cubic(shape=[3,3,3])
>>> geom_1 = OpenPNM.Geometry.Stick_and_Ball(network=pn,pores=[3,6,9],name='geom_1')
>>> geom_1['pore.map']
array([3, 6, 9])

This allows Geometry objects to use to topological queries of the Network class, so that it can find neighbors of its own pore 1, since it knows that maps onto pore 6 of the Network.  

Another important step is the pore and throat labels arrays are created in the Network corresponding to the Geometry name.  So if Geometry called 'geom_1' is associated with a Network, the 'pore.geom_1' and 'throat.geom_1' label arrays are added to Network:

>>> pn.labels()
['pore.all', 'pore.back', 'pore.bottom', 'pore.front', 'pore.geom_1', 'pore.internal', 'pore.left', 'pore.right', 'pore.top', 'throat.all', 'throat.geom_1']
>>> pn.pores('geom_1')
array([3, 6, 9], dtype=int64)

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Adding Models
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
A Geometry object is associated with the network, ``pn``, and applied to all the specified pores and throats in the network.  At this point, however, it is still an empty object with no pore or throat size *information*. To begin assigning size information, it is possible to simply assign values:

>>> geom['pore.seed'] = sp.rand(geom.Np)
>>> geom['throat.constant'] = 1.4

There are, however, very few cases where such a simple assignment is sufficient and usually more elaborate pore scale models will be invoked.  For instance, it is common for throats to adopt the smaller of the seed values in it's two neighboring pores.  OpenPNM includes a library of pre-written models.  Pore scale geometry models are located under ``OpenPNM.Geometry.models``.  There are numerous files in this library with names that indicate their contents (i.e. pore_volume), and each of these files contain a variety of functions for calculating that property.  Specifying which models to use for a given property is done using the ``add_model`` method:

>>> import OpenPNM.Geometry.models as gm
>>> geom.add_model(propname='pore.seed',model=gm.pore_misc.random)

The above call to ``add_model`` does several things.  Firstly, it adds an array to the ``geom`` dictionary called 'pore.seed'.  Secondly, it runs the function it received for the model argument and stores the returned values in 'pore.seed'.  Finally, it saves the model in a private dictionary on the `geom` object.  This final step is essential so that the *Geometry* object can retain a memory of it's models.  This means that the properties of the *Geometry* can be *regenerated* (using ``regenerate``).  This also has the added benefit that the object can be saved to disk and still function fully when it's reloaded. 
	
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Multiple Geometries on a Single Network
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
There can be multiple *Geometry* objects defined for different locations in a *Network* simultaneously.  This was intended to allow for multi-layer media (such fuel cell gas diffusion layers with micro-porous layers on one side), but is also quite useful when applying boundary pores which usually need to have special pore geometry such as 0 volume to produce consistent results.  For instance:

>>> pn = OpenPNM.Network.Cubic(shape=[3,3,3])
>>> Ps1 = pn.pores('top')
>>> geom1 = OpenPNM.Geometry.GenericGeometry(network=pn,pores=Ps1)
>>> Ps2 = pn.pores('top',mode='not')
>>> Ts2 = pn.throats('all')
>>> geom2 = OpenPNM.Geometry.GenericGeometry(network=pn,pores=Ps2,throats=Ts2)

.. note:: **Accessing Geometry Data Via the Network**

    One of the complications that arises from allowing multiple Geometry objects is that the pore size data for the Network becomes distributed across several objects.  This makes it challenging for algorithms to operate on the entire network at once.  To circumvent this problem, the Network object has the special ability to gather Geometry data from all of it's Geometry objects and return them as a single array:

    >>> geom1['pore.seed'] = 0.2
    >>> geom2['pore.seed'] = 0.8
    >>> pn['pore.seed']
    array([ 0.8,  0.8,  0.2,  0.8,  0.8,  0.2,  0.8,  0.8,  0.2,  0.8,  0.8,
            0.2,  0.8,  0.8,  0.2,  0.8,  0.8,  0.2,  0.8,  0.8,  0.2,  0.8,
            0.8,  0.2,  0.8,  0.8,  0.2])

    If any of the Geometry object do not have the requested property, then NaN values are inserted into it's pore/throat locations.  
    
    This special ability is not reversible, meaning that it is not possible to *write* to all Geometry objects from Network:
    
    >>> pn['pore.seed'] = 0.5

    Attempting to do so will result in the error "pore.seed is already defined in at least one associated Geometry object".
		
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Customizing Geometry
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For description of how to create customized subclasses, add properties to the model library, and add new models see :ref:`Customizing OpenPNM<customizing>`











