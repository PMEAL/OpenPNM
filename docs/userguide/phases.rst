.. _phases:

###############################################################################
Phases
###############################################################################
The *Phase* module controls the physical properties of the phases used inside the pore network.  Fluid properties such as gas viscosity and liquid density are calculate by the *Phase* classes.  These are necessary when performing quantitative simulations of transport processes in the network.  The *Phase* module works very similar to the *Geometry* model outlined above.  There is a ``GenericPhase`` class that is associated with a network object during instantiation as follows:

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

.. note:: Phases Exist Everywhere
	
	Notice that pores and throats were *not* sent to the GenericPhase constructor.  This is because *Phases* exist everywhere.  This might seem counterintuitive in a multiphase simulation where one phase displaces another, but it is much easier to calculate the *Phase* properties everywhere, and separately track where each phase is present and in what amount.  

The pore temperature and pressure are automatically set to standard conditions when the object is generated.  These can easily be changed:

>>> air['pore.temperature'] = 363.0

Note that sending a scalar (363.0) without including an index into 'pore.temperature' indicates that all pore temperatures should be set to the same value.  Specific pore values can also be set with:

>>> air['pore.temperature'][[0,1,2,3]] = 355.0

Phase objects are meant to do more than just track temperature.  Many models for estimating the various thermo-physical properties of fluids and solids have been included.  These are located in OpenPNM.Phases.models, with each file corresponding to a common property, containing several options for gas, liquid and solid properties.  Models are retrieved from this library and attached to the *Phase* objects as follows:

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

**Phase** objects are 'built' by the user to contain the specific methods that are to be used to calculate the phase properties.  For instance, a **Phase** object can calculate viscosity assuming a constant value, or Reynolds equation.  The user can use the methods supplied with OpenPNM, or add their own.  

Typically there will be multiple **Phase** objects defined for each simulation, since most models will have at least an invading fluid and a defending fluid.  There can be an unlimited number of phases associated with a **Network**.  
