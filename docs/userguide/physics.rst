.. _physics:

===============================================================================
Physics
===============================================================================
Physics objects are where geometric data and fluid property data are combined to compute the pore-scale physical behavior in the simulation.  For instance, the capillary entry pressure for a throat is a function of size (from Geometry) and the surface tension of the fluid (from Phase), but there are many ways to compute the actual entry pressure, including the Washburn equation for a cylinder, or the Purcell equation for a toroid.  Specifying unique pore-scale Physics models is what sets pore network simulations apart from each other.  The Physics object manages these pore-scale properties and models.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Basic Usage
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Physics objects have the most requirements for initialization.  A Physics object must be associated with a Network and a Phase, and must also be assigned to certain pores and throats.

.. code-block:: python

	pn = OpenPNM.Network.Cubic(shape=[3,3,3])
	Ps = pn.pores('all')
	Ts = pn.throats('all')
	air = OpenPNM.Phases.Air(network=pn)
	phys = OpenPNM.Physics.Standard(network=pn,phase=air,pores=Ps,throats=Ts,name='phys_1')
	
Several important event occur upon instantiation of a Physics object.  Firstly, 'pore.map' and 'throat.map' properties are created in the Physics dictionary.  These contains the pore and throat number in the main Network where the Physics is defined.  Secondly, label arrays are created in the Phase object, called 'pore.phys_1' and 'throat.phys_1'.  These labels indicate where in the full Network the Physics object applies.  Between the 'map' arrays and 'label' arrays it is simple to translate the pore and throat numbers between Physics to the main Network.  

It is possible to set the pore and throat locations of a Physics object after instantiation.  The ``set_locations`` method accepts a list of pores and/or throats and then updates the 'map' and label arrays accordingly.  This only works on an object that has NOT previously been assigned to any pores and/or throats.  Once the assignment is made it cannot be undone.  This functionality is useful for associating a saved Physics objects with a simulation.  

-------------------------------------------------------------------------------
Accessing Physics Data Via the Phase
-------------------------------------------------------------------------------
Each Physics object must be associated with a single Phase object since a physical property like hydraulic conductance depends on the Phase viscosity, which is unique to a Phase.  Physics objects, however, do NOT need to be associated with a Geometry object despite the fact that a physical property like hydraulic conductance also depends on throat size.  The Physics objects use the ability of Network to query all the Geometry objects to extract the requested information from the specified pores and throats, even if they span across multiple Geometries. 

Like Geometry object, Physics objects can lead to a sort of 'fragmentation' of data across multiple objects.  It is ofter desired to retrieve *all* the Physics properties for the entire Network in a single call.  For instance, in fluid flow simulations the hydraulic conductance of the entire Network is needed to construct the coefficient matrix in the solver.  To remedy this, Phase objects have the special ability to query all of their respective Physics objects and retrieve properties from any or all pore and throat locations.  For instance:

>>> phys1 = OpenPNM.Physics.GenericPhysics(network=pn,phase=air,pores=pn.pores('top'))
>>> phys2 = OpenPNM.Physics.GenericPhysics(network=pn,phase=air,pores=pn.pores('top',mode='not'))
>>> phys1['pore.test_vals'] = 1.0
>>> phys2['pore.test_vals'] = 2.0
>>> air['pore.test_vals']
array([ 2.,  2.,  1.,  2.,  2.,  1.,  2.,  2.,  1.,  2.,  2.,  1.,  2., 2.,  1.,  2.,  2.,  1.,  2.,  2.,  1.,  2.,  2.,  1.,  2.,  2.,  1.])

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Customizing Physics
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For description of how to create customized subclasses, add properties to the model library, and add new models see :ref:`Customizing OpenPNM<customizing>`