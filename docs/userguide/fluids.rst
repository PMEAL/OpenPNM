.. _fluids:

###############################################################################
Fluids and Phases Module
###############################################################################
Virtually all algorithms in OpenPNM **require** fluid properties of some sort, either directly or through the pore scale physics methods.  The also **produce** new fluid properties.  The various roles of this object are outlined below.

===============================================================================
What is a Fluid in OpenPNM?
===============================================================================
A fluid offers 4 functions:

1. All the physical properties of the fluid, such as viscosity.  Each property is stored in a dictionary whose key is the property name.  For instance, the viscosity of water in the pores would be accessed with `water.pore_conditions['viscosity']`.  Similarly, throat conditions, if they exist are stored as `water.throat_conditions['viscosity']`.  If the property is uniform throughout the network, then the value can be stored as a scalar to save memory and OpenPNM will apply it correctly.  

2. In a manner similar to first function, the Fluid object also stores the results of and Algorithms and pore-scale Physics calculations that are performed.  The fluid object thus embodies the outcome of the series of simulations that have been run with it.  

3. Functions for altering the fluid state, such as refreshing or resetting conditions.

4. A dictionary containing the *recipe* that was used to define the fluid.  The creation of a fluid object requires a dictionary that contains the necessary parameters which include some key physical characteristics of the fluid, like critical temperature, the names of the property estimation methods desired, and the necessary parameters that should be used in the desired methods.  

-------------------------------------------------------------------------------
Creating a Fluid
-------------------------------------------------------------------------------
The parameters dictionary that must be given to the `create()` looks like:

.. code-block:: python

	water_recipe = {     
		'name': 'water',
		'Pc': 2.206e6, #Pa
		'Tc': 647,     #K
		'MW': 0.0181,  #kg/mol
		'diffusivity':{
			'method': 'constant',
			'value': 1e-12},
		'viscosity':{
			'method': 'constant',
			'value': 0.001},
		'molar_density':{
			'method': 'constant',
			'value': 44445},
		'surface_tension':{
			'method': 'Eotvos',
			'k': 2.25e-4},
		'contact_angle':{
			'method': 'constant',
			'value': 120},
	}

The above example used the `constant` method for all properties, but many other methods are available as listed below.


-------------------------------------------------------------------------------
Refreshing a Fluid for Changes in Network Conditions
-------------------------------------------------------------------------------
.. automethod:: OpenPNM.Fluids.GenericFluid.regenerate

For example, refreshing water properties for changes in the network can be achieved as follows:

.. code-block:: python

	water.regenerate()

-------------------------------------------------------------------------------
Available Property Estimation Methods
-------------------------------------------------------------------------------
.. automodule:: OpenPNM.Fluids.Diffusivity
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Fluids.Viscosity
   :members:
   :undoc-members:
   :show-inheritance:
 
.. automodule:: OpenPNM.Fluids.MolarDensity
   :members:
   :undoc-members:
   :show-inheritance:
   
.. automodule:: OpenPNM.Fluids.SurfaceTension
   :members:
   :undoc-members:
   :show-inheritance:
   
.. automodule:: OpenPNM.Fluids.ContactAngle
   :members:
   :undoc-members:
   :show-inheritance:
 
.. automodule:: OpenPNM.Fluids.ThermalConductivity
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Fluids.VaporPressure
   :members:
   :undoc-members:
   :show-inheritance:

