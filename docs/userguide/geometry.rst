.. _geometry:

###############################################################################
Geometry
###############################################################################
The *Geometry* module manages the network pore and throat size information.  This module contains the ``GenericGeometry`` class, which like all OpenPNM objects is subclass of Python's ``dict`` class, but has numerous OpenPNM specific methods added to it.  An empty ``GenericGeometry`` object *can* be initialized with no arguments, but this is not a useful object since it isn't associated with a network or assigned to any pores.  A more useful *Geometry* object is obtained by instantiating a non-empty network, then assign a GenericGeometry to all pores and throats:

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

This Geometry object is now associated with the network, ``pn``, and applied to all the pores and throats in the network.  At this point, however, it is still an empty object with no pore or throat size *information*.  To begin assigning size information, it is possible to simply assign values:

>>> geom['pore.seed'] = sp.rand(geom.Np)
>>> geom['throat.constant'] = 1.4

There are, however, very few cases where such a simple assignment is sufficient and usually more elaborate pore scale models will be invoked.  For instance, it is common for throats to adopt the smaller of the seed values in it's two neighboring pores.  OpenPNM includes a library of pre-written models.  Pore scale geometry models are located under ``OpenPNM.Geometry.models``.  There are numerous files in this library with names that indicate their contents (i.e. pore_volume), and each of these files contain a variety of functions for calculating that property.  Specifying which models to use for a given property is done using the ``add_model`` method:

>>> import OpenPNM.Geometry.models as gm
>>> geom.add_model(propname='pore.seed',model=gm.pore_misc.random)

The above call to ``add_model`` does several things.  Firstly, it adds an array to the ``geom`` dictionary called 'pore.seed'.  Secondly, it runs the function it received for the model argument and stores the returned values in 'pore.seed'.  Finally, it saves the model in a private dictionary on the `geom` object.  This final step is essential so that the *Geometry* object can retain a memory of it's models.  This means that the properties of the *Geometry* can be *regenerated* (using ``regenerate``).  This also has the added benefit that the object can be saved to disk and still function fully when it's reloaded.

.. note:: Object Inheritance
	
	Because all objects in OpenPNM descend from the same base class, they all have many of the same query methods.  You can query the number of pores (``num_pores``) and throats (``num_throats``), find a list of available labels (``labels``) and properties (``props``), determine the name of the object (``name``), and retrieve a list of all the other objects associated with the same network (``geometries``, ``physics``, ``phases``).  

There can be multiple *Geometry* objects defined for different locations in a *Network* simultaneously.  This was intended to allow for multi-layer media (such fuel cell gas diffusion layers with micro-porous layers on one side), but is also quite useful when applying boundary pores which usually need to have special pore geometry such as 0 volume to produce consistent results.  For instance:

>>> pn = OpenPNM.Network.Cubic(shape=[3,3,3])
>>> Ps1 = pn.pores('top')
>>> geom1 = OpenPNM.Geometry.GenericNetwork(network=pn,pores=Ps1)
>>> Ps2 = pn.pores('top',mode='not')
>>> Ts2 = pn.throats('all')
>>> geom2 = OpenPNM.Geometry.GenericNetwork(network=pn,pores=Ps2,throats=Ts2)

It is of possible to develop customized pore-scale models beyond those included with OpenPNM.  This is outlined in more detail :ref:`here<geometry>`.
