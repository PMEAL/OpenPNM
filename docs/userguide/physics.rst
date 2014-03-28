.. _physics:

###############################################################################
Pore Scale Physics
###############################################################################

.. note:: 

	Fluid, Geometry and Physics modules are designed to function identically, so once you're familiar with the usage of one then all the others should be similar.  
	
===============================================================================
What is a Physics Object?
===============================================================================

===============================================================================
Generating Physics Objects
===============================================================================
blah

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Building a Physics Object
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
blah

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Generating or Regenerating Physics Data
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

===============================================================================
Customizing the Physics Submodules
===============================================================================

.. _custom_prop_names:
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Handling Custom Fluid and Geometry Property Names
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
When a Physics method is run, it assumes the requisite properties are stored in the default locations.  For instance, it assumes that the diffusivity of a fluid can be found with fluid.get_pore_data(prop='diffusivity'). If, however, a custom name has been used when adding the 'diffusivity' property to the fluid, then the data will be stored with that custom property name.  To handle this eventuality, all Physics methods accept the property name as optional arguments.  For example, calculating the diffusive conductance of fluid 1 using a custom property name of 'DAB' would be done as:

.. code-block:: python

	phys.add_method(prop='diffusive_conductance',model='fuller',diffusivity='DAB',**kwargs)

Where ``**kwargs`` represents the additional parameters required by the fuller model.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Adding Custom Models to Existing Physics Classes
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This documentation is being rewritten, sorry for the inconvenience.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Adding Custom Models to Existing Physics Classes
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This documentation is being rewritten, sorry for the inconvenience.


===============================================================================
Available Property Estimation Models
===============================================================================

For a complete list of available pore scale physics models see the :ref:`Function Reference <physics_ref>`.