.. _phases:

===============================================================================
Phases
===============================================================================
The *Phase* module controls the physical properties of the phases used in simulation.  Thermophysical properties such as gas viscosity and liquid density are calculated by the *Phase* objects.  These are necessary when performing quantitative simulations of transport processes in the network.  

.. inheritance-diagram:: OpenPNM.Phases.GenericPhase

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Basic Usage
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
There is a ``GenericPhase`` class that is associated with a Network object during instantiation as follows:

>>> air = OpenPNM.Phases.GenericPhase(network=pn,name='air')
>>> print(air)
------------------------------------------------------------
OpenPNM.Phases.GenericPhase: 	air
------------------------------------------------------------
#     Properties                          Valid Values
------------------------------------------------------------
1     pore.pressure                          27 / 27   
2     pore.temperature                       27 / 27   
------------------------------------------------------------
#     Labels                              Assigned Locations
------------------------------------------------------------
1     pore.all                            27        
2     throat.all                          54        
------------------------------------------------------------

Note that the Phase is defined everywhere in the Network:

>>> pn.num_pores()
27
>>> air.num_pores()
27

Unlike Geometry objects, Phases do not need to create any 'map' arrays or label arrays in the Network since Phases are applied everywhere.  Although it might seem counter-intuitive in a multiphase simulation to have Phases existing everywhere, it is simply a matter of convenience to calculate the Phase properties for the entire domain.  Tracking the actual locations of the Phases is handled separately, usually by an Algorithm object that performs some sort of percolation calculation.  

Typically there will be multiple **Phase** objects defined for each simulation, since most models will have at least an invading fluid and a defending fluid.  There can be an unlimited number of phases associated with a **Network**.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Adding Models
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The pore temperature and pressure are automatically set to standard conditions when the object is generated.  These can easily be changed:

>>> air['pore.temperature'] = 363.0

Note that sending a scalar (363.0) without including an index into 'pore.temperature' indicates that all pore temperatures in the Network should be set to the same value.  Specific pore values can also be set with:

>>> air['pore.temperature'][[0,1,2,3]] = 355.0

Phase objects are meant to do more than just track temperature.  Many models for estimating various thermo-physical properties of fluids and solids have been included.  These are located in OpenPNM.Phases.models, with each file corresponding to a common property, containing several options for gas, liquid and solid properties.  Models are retrieved from this library and attached to the *Phase* objects as follows:

>>> air['pore.temperature'] = 298.0
>>> model = OpenPNM.Phases.models.molar_density.ideal_gas
>>> air.add_model(propname='pore.molar_density',model=model)
>>> air['pore.molar_density'][0]
40.894621255593904

Now if the temperature of the pores is changed, all the other properties will also be changed, provided that the other property models are coded as functions of temperature:

>>> air['pore.temperature'] = 353.0
>>> air.regenerate()
>>> air['pore.molar_density'][0]
34.522938057130261

Note that the ``regenerate`` method must be called for the change in temperature to propagate to the other properties.  

.. note:: The Meaning of Pore vs Throat Properties in Phases

    In general all Phase properties are specified in pores only.  The basis of this is that most algorithms solve for the conditions in the pores.  For instance, a FourierConduction algorithm solves for pore temperatures, so a temperature dependent viscosity should also be evaluated in the pores.  This holds for most properties.  
	
	Ironically, in order to solve for properties in the pores it is usually necessary to know the Phase conditions in the throats.  OpenPNM objects include the ability to interpolate throat conditions based on the conditions in the neighboring pores.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Mixtures
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
In many cases of practical interest the 'phase' is actually a mixture.  It is possible to create a mixture Phase by simply sending the 'pure component' phases as arguments to the initialization of the mixture:

.. code-block:: python

    N2 = OpenPNM.Phases.GenericPhase(network=pn,name='pure_N2')
    O2 = OpenPNM.Phases.GenericPhase(network=pn,name='pure_O2')
    air = OpenPNM.Phases.GenericPhase(network=pn,name='air',components=[N2,O2])

The key difference between the instantiation of 'N2' and 'O2' versus 'air' is that 'air' receives the others as 'components'.  During the initialization of a Phase any Phases received as 'components' are associated with the mixture as can be seen with:

>>> air.phases()
['pure_O2','pure_N2']

With this association it is now possible to extract pure component property information from each component Phase which can be used to calculate mixture properties.  There is one caveat with this approach however: the 'composition' of each component in the mixture (i.e. mole fraction of each component) is stored on the individual component Phases under the 'pore.mole_fraction' property name.  The mole fraction of each Phase can be specified as:

.. code-block:: python

    N2['pore.mole_fraction'] = 0.79
	O2['pore.mole_fraction'] = 0.21
	N2['pore.molar_mass'] = 0.028
	O2['pore.molar_mass'] = 0.032

With this information it is possible to calculate mixture properties such as the average molecular weight and so on.  There are a small number of Phase models to work with mixtures at present, found throughout the various fluid property model categories:

.. code-block:: python
	
    mod = OpenPNM.Phases.models.molar_mass.mixture
    air.add_model(propname='pore.molar_mass',model=mod)
    air['pore.molar_mass'][0]
    0.02884
	
The ``mixture`` molar mass method looks into the Phases of the mixture, retrieves their molar masses and mole fractions, and computes the average molar mass of the mixture.  	

.. note:: Mixture Temperature and Pressure

    Temperature and pressure are the two thermodynamic properties required for calculating most other phase properties.  As the mixture temperature and pressure change, it is necessary to also update the temperature and pressure of the component phases so their properties are calculated correctly.  When a component Phase is associated with a mixture, OpenPNM automatically adds a model to the component Phase that forces its temperature and pressure to match that of the mixture.  Thus, whenever the mixture temperature or pressure change, so do the component Phases.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Customizing Phases
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For description of how to create customized subclasses, add properties to the model library, and add new models see :ref:`Customizing OpenPNM<customizing>`
