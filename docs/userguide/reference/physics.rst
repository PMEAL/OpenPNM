.. _physics:

===============================================================================
Physics
===============================================================================
**Physics** objects are where geometric data and phase property data are combined to compute the pore-scale physical behavior within the simulation.  For instance, the capillary entry pressure for a throat is a function of size (from a **Geometry**) and the surface tension of the fluid (from a **Phase**), but there are many ways to compute the actual entry pressure, including the *Washburn equation* for a cylinder, or the *Purcell equation* for a toroid.  Specifying unique pore-scale Physics models is what sets pore network simulations apart from each other.  The Physics object manages these pore-scale properties.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Basic Usage
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
**Physics** objects have the most requirements for initialization.  A **Physics** object must be associated with a **Network**, a **Geometry** object that implies which pores and throats the **Physics** applies to, and a **Phase** object which specifies which thermophysical properties apply.

.. code-block:: python

	pn = OpenPNM.Network.Cubic(shape=[3,3,3])
	Ps = pn.pores('all')
	Ts = pn.throats('all')
	geom = OpenPNM.Geometry.Stick_and_Ball(network=pn,pores=Ps,throats=Ts)
	air = OpenPNM.Phases.Air(network=pn)
	phys = OpenPNM.Physics.Standard(network=pn,phase=air,pores=Ps,throats=Ts,name='phys_1')

Several important events occur upon instantiation of a Physics object.  Mainly, *label* arrays are created in the **Phase** object, called 'pore.phys_1' and 'throat.phys_1'.  These indicate which pores and throats in full domain the **Physics** object applies.  These locations are recorded in the **Phase** object to enable data exchange between the **Phase** and **Physics**, as described in the following note.

.. note:: **Accessing Physics Data Via the Phase**

	Each **Physics** object must be associated with a single **Phase** object since a physical property like hydraulic conductance depends on the **Phase** viscosity, which is unique to a **Phase**.  
	
	Like a **Geometry** object, the encapsulation of data on individual **Physics** objects can lead to a sort of *fragmentation* of data across multiple objects.  It is often desired to retrieve *all* the **Physics** properties for the entire **Network** in a single call.  For instance, in fluid flow simulations the hydraulic conductance of the entire Network is needed to construct the coefficient matrix in the solver.  To remedy this, **Phase** objects have the special ability to query all of their respective **Physics** objects and retrieve properties from any or all pore and throat locations.  For instance:
	
	>>> phys1 = OpenPNM.Physics.GenericPhysics(network=pn,phase=air,pores=pn.pores('top'))
	>>> phys2 = OpenPNM.Physics.GenericPhysics(network=pn,phase=air,pores=pn.pores('top',mode='not'))
	>>> phys1['pore.test_vals'] = 1.0
	>>> phys2['pore.test_vals'] = 2.0
	>>> air['pore.test_vals']
	array([ 2.,  2.,  1.,  2.,  2.,  1.,  2.,  2.,  1.,  2.,  2.,  1.,  2., 2.,  1.,  2.,  2.,  1.,  2.,  2.,  1.,  2.,  2.,  1.,  2.,  2.,  1.])

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Adding Models
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
OpenPNM includes a collection of prewritten pore-scale physics models.  They are located under *'Physics.models'*.  They are further categorized according to the phenomena they describe (i.e. *'Physics.models.capillary_pressure'*).  Inside each file there are several options (e.g. *washburn*, and *purcell*).  To utilize one of these models on a **Physics** object it must be added using ``add_model``:

>>> model = OpenPNM.Physics.models.capillary_pressure.washburn
>>> phys.add_model(propname='throat.capillary_pressure',model=model)

This command tells ``phys`` that you want to use the *'washburn'* model to calculate throat entry pressures, and that you would like the values returned by this function to be stored under the property name *'throat.capillary_pressure'*.  

Each Physics model assumes the pre-existence of the necessary **Phase** property data.  The *'capillary_pressure'* models require that the *'surface tension'* and *'contact angle'* are present on the **Phase** object that the **Physics** is associated with.  To learn about using non-default property names, see :ref:`Customizing OpenPNM<customizing>`

When ``add_model`` is called, the model being added is run and the data are stored in the Physics dictionary under the specified *'propname'*.  If some conditions in the simulation change, like *'temperature'*, then this will impact the **Phase** properties (i.e. *'throat.surface_tension'*).  For such changes to propagate through the simulation the ``regenerate`` method must be called.  If the **Phase** temperature did change, then ``Phase.regenerate`` will update the *'surface_tension'* values (assuming that a temperature dependent model was assigned to calculate it).  It is then necessary to call ``Physics.regenerate`` so that the *'capillary_pressure'* values are recalculated to reflect the changed *'surface_tension'* values.  

For a full description of adding and manipulating models see :ref:`Pore Scale Models<models>`.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Customizing Physics
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For description of how to create customized subclasses, add properties to the model library, and add new models see :ref:`Customizing OpenPNM<customizing>`