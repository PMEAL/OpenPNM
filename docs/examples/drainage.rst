.. _drainage-example:

===============================================================================
Drainage Curve on a Cubic Network
===============================================================================

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Generating the Network, adding Geometry and creating Phases
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Start by generating a basic cubic network and the other required components:

.. code-block:: python

	import OpenPNM
	pn = OpenPNM.Network.Cubic(name='net',shape=[10,10,10],spacing=0.0001)
	pn.add_boundaries()

The last call adds a layer of boundary pores around the network after it is generated. These boundary pores will be used in the following calculations. Next we generate a geometry for the network and the phase, in this case air. A geometry can span over a part of the network only, so we need to specify to which pores and throats this geometry object should apply. For this example, we want it to apply to all pores and throats of the network. To do so, we can get all pore and throat indices of the network with the ``pn.pores()`` and ``pn.throats()`` calls, and pass these to the geometry object.

.. code-block:: python

	Ps = pn.pores()
	Ts = pn.throats()
	geom = OpenPNM.Geometry.Toray090(network=pn,name='basic',pores=Ps,throats=Ts)
	air = OpenPNM.Phases.Air(network=pn)
	water = OpenPNM.Phases.Water(network=pn)

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Define the Pore-scale Physics
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

To perform most algorithms, it is necessary to define the pore scale physics that relates pore/throat geometry, phase properties, and the mechanism to be modeled.  In the case of a drainage curve simulation, it is necessary to define the pressure at which the invading phase can enter each throat.  This is commonly done by assuming the throat is a cylinder and using the so-called 'Washburn' equation.  The OpenPNM Physics module has a submodule for capillary pressure methods, including the Washburn model.  To use this model in a simulation, you first create a generic Physics object. 

>>>	phys = OpenPNM.Physics.GenericPhysics(network=pn,phase=water,geometry=geom)

Then add the desired methods to this object using:

>>>	phys.add_model(propname='throat.capillary_pressure',model=OpenPNM.Physics.models.capillary_pressure.washburn)

This means that the Physics object will now have a function called 'capillary_pressure', that when called will calculate throat entry pressures using the 'washburn' model.  The Washburn model requires that the Phase object (Water in this case) has the necessary physical properties of surface tension and contact angle.

.. note::

	Both surface tension and contact angle are actually 'phase system' properties, rather than solely water properties.  It is an open problem in OpenPNM to figure out how to treat these sort of properties more rigorously.  For the present time, they must be entered a single phase properties.

The predefined Water object is assigned a contact angle of 110 degrees by default (water on Teflon). This can be confirmed by first printing a list of which properties have been assigned to Water, and then checking the values.

>>> print(water)
------------------------------------------------------------
1 	:  pore.contact_angle
2 	:  pore.critical_pressure
3 	:  pore.critical_temperature
4 	:  pore.critical_volume
5 	:  pore.density
6 	:  pore.diffusivity
7 	:  pore.molar_density
8 	:  pore.molecular_weight
9 	:  pore.pressure
10 	:  pore.surface_tension
11 	:  pore.temperature
12 	:  pore.thermal_conductivity
13 	:  pore.vapor_pressure
14 	:  pore.viscosity
------------------------------------------------------------
>>> print(water['pore.contact_angle'])[0]
110.0

To change this value, the 'pore.contact_angle' property of the Water object can be set to a new constant:

.. code-block:: python

    water['pore.contact_angle'] = 140.0

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Run a Drainage Simulation
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

At this point, the system is fully defined and ready to perform some simulations.  A typical algorithm used in pore network modeling is to use ordinary percolation to simulate drainage of wetting phase by invasion of a nonwetting phase.  An Algorithm object is be created as follows:

>>>	OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=pn,invading_phase=water)

Before performing simulations with this algorithm it is necessary to specify the desired experimental parameters in the ``run()`` command:

>>>	Ps = pn.pores(labels=['bottom_boundary'])
>>>	OP_1.run(inlets=Ps)

The first line in the finds all the pores in the network that are labeled 'bottom_boundary' and assigns it to 'Ps'.  This labeling step was applied during the network construction.  The list of pores which are to be considered as phase inlets along with which phases are the invading and defending phase are set to the `run()` method and the algorithm proceeds.  Upon completion one can view resulting capillary pressure curving using the following command:

>>>	OP_1.plot_drainage_curve()

The red and blue lines represent the filling of pores and throats separately.  The non-zero starting point of the red lines (pores) is due to the fact that the inlet pores are invaded at the start of the process.  This can be avoided by defining a second geometry for boundary pores that have zero volume.  

-------------------------------------------------------------------------------
Sharing Algorithm Results Throughout the Simulation
-------------------------------------------------------------------------------

The results of the above simulation (and all simulations) are stored locally on the algorithm object.  If these results are to be used in other parts of the simulations, then they must be explicitly sent 'out'.  Keeping the results *silo-ed* in this way prevents unintentional overwriting of results by subsequent algorithms.  This allows for multiple simulations of the same type to be run with different conditions and such.  Sending the results of any simulation 'out' is done by with the `return_results()` command:

>>>	OP_1.return_results(Pc=8000)

The above command outputs data called 'occupancy' to the invading phase object. This data describes which pores and throats are filled by invading and defending phase at the specified applied capillary pressure *Pc*.  This information can be used by subsequent algorithms.  For instance it is often of interest to determine the gas phase diffusivity through a partially water filled network.  The Fickian diffusion algorithm then would use this information and set gas diffusion through water filled pores to zero and a relative effective diffusivity value could be found.