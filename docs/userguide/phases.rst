.. _phases:

===============================================================================
Phases
===============================================================================
The *Phase* module controls the physical properties of the phases used in simulation.  Thermophysical properties such as gas viscosity and liquid density are calculated by the *Phase* objects.  These are necessary when performing quantitative simulations of transport processes in the network.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Basic Usage
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
There is a ``GenericPhase`` class that is associated with a Network object during instantiation as follows:

>>> air = OpenPNM.Phases.GenericPhase(network=pn,name='air')
>>> print(air)
------------------------------------------------------------
OpenPNM.Phases.GenericPhase: 	GenericPhase_jypSa
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
54

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

Note that the ``regenerate`` method must called for the change in temperature to propagate to the other properties.  

Typically there will be multiple **Phase** objects defined for each simulation, since most models will have at least an invading fluid and a defending fluid.  There can be an unlimited number of phases associated with a **Network**.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Customizing Phases
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For description of how to create customized subclasses, add properties to the model library, and add new models see :ref:`Customizing OpenPNM<customizing>`
