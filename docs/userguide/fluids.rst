.. _fluids:

###############################################################################
Fluid Properties
###############################################################################
Virtually all algorithms in OpenPNM require **fluid** properties of some sort, either directly or through the pore scale physics methods.  Many algorithm also produce new fluid properties.  Generating pure fluid properties, and storing the results of pore scale physics and algorithms is performed by the Fluids object, as outlined below.  

.. note:: 

Fluid, Geometry and Physics modules are designed to function identically, so once you're familiar with the usage of one then all the others should be similar.  

===============================================================================
What is a Fluid Object in OpenPNM?
===============================================================================
Fluid objects have 2 main functions in OpenPNM:

1. They possess methods for calculating physical properties of fluids as a function of relevant conditions.  For instance, the Fluid object may possess an attribute called `viscosity`, which when called calculates the viscosity of the fluid based on the temperature of the network.  If there happens to be a temperature distribution in the network, then this method will calculate the resulting distribution of viscosity.  

2. They store all data pertaining to the fluid.  For instance, the viscosity of the fluid will be stored on the Fluid object.  Additionally, any data indirectly pertaining to fluid behavior is also stored on the fluid, such as hydraulic conductance, which is a product of both the fluid viscosity, but also the pore and throat dimensions.  This information is stored on the fluid in order to compartmentalize data pertaining to each fluid, since gas and liquid will both have hydraulic conductances.  

===============================================================================
Creating a Fluid
===============================================================================
Fluid objects are designed to be highly customizable.  The general process of creating a fluid involves first initializing the Fluid object as shown below.  Note that the initialization takes a pore network object as an argument.  This is necessary so the Fluid is aware of the network topology and size, so it can store data for each pore and throat.  

.. code-block:: python

  >>> pn = OpenPNM.Network.TestNet()  # Creates a simple 5 x 5 x 5 network for testing
  >>> liq_1 = OpenPNM.Fluids.GenericFluid(network=pn,name='liquid')
  
This basic object does not contain any property estimation methods.  These must be selected and added individually.  OpenPNM includes a number of submodules under the Fluids module, such as `viscosity`, `diffusivity`, `molar_density` and so on.  Each of these submodules has multiple predefined models available for calculating each property.  For instance the `molar_density` submodule has a model called `ideal_gas_law`.  Most or all of these models take input arguments that customize their output to a specific fluid.  The `ideal_gas_law` model requires only R (the gas constant in the appropriate units); the `viscosity` submodule, on the other hand, has the `reynolds` model which requires uo and b that are fluid specific.  Methods are added to the Fluid object using `add_method` as follows

.. code-block:: python

  >>> liq_1.add_method(prop='molar_density', model='constant', value=100)  # apply a constant density
  >>> liq_1.add_method(prop='viscosity', model='reynolds', uo=1, b=1)  # use a temperature dependent model
  >>> liq_1.add_method(prop='molar_mass', prop_name='MW', model='constant', value=0.02)  # use a custom property name

The methods are added to the Fluid object according to the property name (`prop`) by default, but can be given customized names as well using the `prop_name` argument.  The number of fluid properties that are added is arbitrary and customizable.  If only fluid flow calculations will be performed, then it is not necessary to add a `diffusivity` method to the fluid.  Once the desired methods have been added, the next step is to actually calculate the fluid properties.  This is done by calling the added methods as follows:

.. code-block:: python

  >>> liq_1.molar_density()
  >>> liq_1.viscosity()

Whenever the data need to be updated, such as when the temperature of the network changes, then the methods need to be called again.  The Fluid object contains a helper function called `regenerate` which will call all of the added methods in the order they were added.  It is also possible to update only certain methods by sending their names as string arguments to `regenerate`.

.. code-block:: python

  >>> liq_1.regenerate()
  >>> liq_1.regenerate(['molar_density','viscosity'])
  
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Adding Custom Property Models
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
OpenPNM comes with a small set of property models for each of the property submodules.  It was designed to be as simple as possible to add new property models to this set.  This is done by simply adding a new method to submodule file of interest.  For example, to add the Peng-Robinson equation of state to the `molar_density` submodule, you simply open OpenPNM/Fluids/molar_density.py and add the function.  There is one caveat: the data produced by the function should be written using the OpenPNM setter method.  This ensures that date is written to the correct location and in the correct format.  It also ensures that the data can be found using the corresponding getter method.  Writing data directly to the Fluid object dictionary is possible, but highly discouraged.  An example of the setter method can be found in any of the provided property model functions.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Adding Custom Properties
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
OpenPNM includes fluid submodules for many common properties, but this list is not exhaustive.  Adding a new property submodule is a two step process.  Firstly, one must create a new file in OpenPNM/Fluids with the desired property name (e.g. enthalpy.py).  Secondly, this file must be added to the __init__.py file in the Fluids folder or else its methods won't be available.  Examples can be found in the __init__.py file, but the required line would be `from . import enthalpy`.  

===============================================================================
Sub-classing a Fluid
===============================================================================
There are several fluids that are used commonly enough that entering individual methods as described above would be annoyingly repetitive.  For these cases, it is helpful to create a subclass of the GenericFluid class that contains a pre-written list of methods to add and the appropriate arguments.  OpenPNM includes subclasses for Air and Water, and these can be used as examples for develop custom subclasses.  There are two steps required to add a custom subclass.  First, a file must be added to the Fluids folder, such as __Oil__.py.  In the initialization method of this file, the various `add_method()` lines that are required to generate the fluid shoudl be added.  Secondly, the new file must be added to the OpenPNM/Fluids/__init__.py file as `from . import __Oil__.py`.  







