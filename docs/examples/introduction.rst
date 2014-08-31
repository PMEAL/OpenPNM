.. _tutorial:

###############################################################################
Introduction: Step by Step Guide to a Standard OpenPNM Script
###############################################################################

===============================================================================
Building a Cubic Network
===============================================================================

The first thing you must do is import the OpenPNM code so you have access to the functions and methods, so in a blank *.py* file or at the python command line, start by entering the following line:

.. code-block:: python
    
    	import OpenPNM
   
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Initialize the Network Topology
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Next, it's time to generate a Network.  This is accomplished by choosing the desired network topology (e.g. cubic), then calling its respective method in OpenPNM with the desired parameters:

.. code-block:: python

	pn = OpenPNM.Network.Cubic(name='net',shape=[10,10,10])

This generates a topological network called *pn* which contains pores at the correct spatial positions and connections between the pores according the desired topology, but without boundary pores.  The network can be queried for certain topological information such as:

.. code-block:: python

	pn.num_pores()  # 1000
	pn.num_throats()  # 2700
	pn.find_neighbor_pores(pores=[1])  # [0,2,11,101]
	pn.labels(pores=[1])  # ['all','bottom','left']
	pn.pores(labels = 'bottom')  # [0,1,2,3,4,5,6,7,8,9]
	pn.throats(labels = 'left')  # [0, 2, 3, 5, 6, .......]

This data may also be stored in a variable:

.. code-block:: python

	Ps = pn.pores()
	Ts = pn.throats()

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Initialize and Build a Geometry Object
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The network does not contain any information about pore and throat sizes at this point.  The next step is to create a geometry object to calculate the desired geometrical properties.  

.. code-block:: python

	geom = OpenPNM.Geometry.GenericGeometry(network=pn,name='stick_and_ball',pores=Ps,throats=Ts)  # instantiate geometry object
	
-------------------------------------------------------------------------------
Add Desired Methods to Geometry
-------------------------------------------------------------------------------
	
This freshly instantiated object contains no methods for actual geometry calculations as yet.  A fully functional object is built by adding the desired methods.  For example, the most basic type of geometry is the so-called 'stick and ball' model, where pores are treated as spheres and throats as cylinders.  Furthermore, it is common to assign pore sizes without regard for spatial correlation, but then to assign throat sizes based on the size of the pores it connects.  This is accomplished by choosing the desired models for each property, then adding them to the geometry object.  

The first step is to load the geometry model library.

.. code-block:: python

	import OpenPNM.Geometry.models as gm

Then, the different geometry models are added one by one to the object geom.

.. code-block:: python

	geom.add_model(propname='pore.seed',model=gm.pore_misc.random,regen_mode = 'static') #begin adding the desired methods to 'geom'
	geom.add_model(propname='throat.seed',model=gm.throat_misc.neighbor,pore_prop='pore.seed',mode='min',regen_mode = 'static')
	geom.add_model(propname='pore.volume',model=gm.pore_volume.sphere,regen_mode = 'static')
	geom.add_model(propname='pore.area',model=gm.pore_area.spherical)
	geom.add_model(propname='throat.length',model=gm.throat_length.straight,regen_mode = 'static')
	geom.add_model(propname='throat.volume',model=gm.throat_volume.cylinder,regen_mode = 'static')
	geom.add_model(propname='throat.area',model=gm.throat_area.cylinder,regen_mode='static')
	
Each of the above commands looks into the submodule associated with the `propname` argument, extracts the model, assigns the specified parameters, and finally attaches the model to the Geometry object.  

OpenPNM ships with many pre-written models available for each property, but adding custom models and even custom properties is designed to be easy.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Create Phases
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

At this point the model is now topologically and geometrically complete.  It has pore coordinates, pore and throat sizes and so on.  In order to perform any simulations, however, it is necessary to build Phases objects that e.g. represent the fluids in the simulations.  This is done using the same composition technique used to build the Geometry.  Phases objects are instantiated and attached to the Network as follows:

.. code-block:: python

	air = OpenPNM.Phases.GenericPhase(network=pn,name='air')
	water = OpenPNM.Phases.GenericPhase(network=pn,name='water')
	
-------------------------------------------------------------------------------
Add Desired Methods to Phases
-------------------------------------------------------------------------------
	
Now it is necessary to fill out these two objects with the desired property calculation model.  For instance, these phases have a very different viscosity and these must be calculated differently.  
As for the geometric object, the phase models need to be load first:

.. code-block:: python

	from OpenPNM.Phases import models as fm

Then, water and air properties are then defined by the code below. Note that some of the models, such as the Fuller model of diffusivity, needs input parameters as molar masses. These inputs are simply state in the add_model method.

.. code-block:: python

	air.add_model(propname='pore.diffusivity',model=fm.diffusivity.fuller,MA=0.03199,MB=0.0291,vA=16.3,vB=19.7)
    	air.add_model(propname='pore.viscosity',model=fm.viscosity.reynolds,uo=0.001,b=0.1)
	air.add_model(propname='pore.molar_density',model=fm.molar_density.ideal_gas,R=8.314)
	water.add_model(propname='pore.diffusivity',model=fm.misc.constant,value=1e-12)
	water.add_model(propname='pore.viscosity',model=fm.misc.constant,value=0.001)
	water.add_model(propname='pore.molar_density',model=fm.misc.constant,value=44445)


	
The above lines retrieve the requested property estimation model from the submodule indicated by the `propname` argument, and assign that method to the corresponding property of the phases on each pore location.  Setting a constant value, as for intance a constant water contact angle, may also be done by directly adding a new dictionnary entry:

.. code-block:: python

	water['pore.contact_angle'] = 110
	water['pore.surface_tension'] = 0.072



+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Create Pore Scale Physics Objects
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

We are still not ready to perform any experiments, despite the fact that phases are defined fully built up.  The last step is to define the desired pore scale physics, which defines how the phase and solid objects interact.  A classic example of this is the Washburn equation which predicts the pressure required to push a non-wetting phase through a capillary of known size.  OpenPNM attempts to permit a high degree of extensibility by using the same object construction approach used for Geometry and Phase above.  Because the Physics object defines the interaction of a Phase with the Geometry, it is necessary to build one physics object for each Phase (and Geometry).  

.. code-block:: python
	phys_water = OpenPNM.Physics.GenericPhysics(network=pn,phase=water,name='standard_water_physics',pores=Ps,throats=Ts)
	phys_air = OpenPNM.Physics.GenericPhysics(network=pn,phase=air,name='standard_air_physics',pores=Ps,throats=Ts)

-------------------------------------------------------------------------------
Add Desired Methods to Physics Objects
-------------------------------------------------------------------------------
	
As with phases and geometry objects, the next steps are first to load the model library and to build-up the bare objects with the desired models:

.. code-block:: python

	from OpenPNM.Physics import models as pm

	phys_water.add_model(propname='throat.capillary_pressure',model=pm.capillary_pressure.purcell,r_toroid=1.e-5)
	phys_water.add_model(propname='throat.hydraulic_conductance',model=pm.hydraulic_conductance.hagen_poiseuille)
	phys_water.add_model(propname='throat.diffusive_conductance',model=pm.diffusive_conductance.bulk_diffusion)

	phys_air.add_model(propname='throat.hydraulic_conductance',model=pm.hydraulic_conductance.hagen_poiseuille) 
	#phys_air.add_model(propname='pore.diffusive_conductance',model='bulk_diffusion')
	phys_air['pore.diffusive_conductance'] = 2e-5
	




+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Visualise the results
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

We can now visualise our geometry and our phase properties. 



-------------------------------------------------------------------------------
Use the Python vtk module
-------------------------------------------------------------------------------

For a quick look, it could be done thanks to the Python vtk module. The following lines below allow you to create the 3D cubic network with spheres 	 representing the pores. The throats are coloured by the value of throats capillary pressure.



.. code-block:: python

	from OpenPNM.Postprocessing.Graphics import Scene, Wires, Spheres
	Cp = water.get_data(prop='capillary_pressure',pores='all',mode='interpolate')
	wires = Wires(pn['pore.coords'], pn['throat.conns'],Cp)
	sphere = Spheres(centers=pn['pore.coords'] ,radii=geom['pore.diameter']*1)  
	scene = Scene()    
	scene.add_actors([wires,sphere])
	scene.play()


-------------------------------------------------------------------------------
Use Paraview
-------------------------------------------------------------------------------
For more detailed visualisaton, the data created by OpenPNM may be exported to a vtk ASCII file to be loaded through Paraview.

.. code-block:: python

	import OpenPNM.Postprocessing.Export as save
	save.VTK(network=pn,phases=[air,water])
	
This creates a *net.vtp* file in the active directory, which can be loaded from ParaView. Visualisation of the pores can be achieved by using 3D Glyphs.