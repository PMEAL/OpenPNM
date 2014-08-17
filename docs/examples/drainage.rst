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
	pn = OpenPNM.Network.Cubic.empty(name='net',loglevel=20,dims=[10,10,10])
	pn.add_boundaries()

In the last call, the property add_boundaries was set to True, which means that a layer of boundary pores around the network is generated. These boundary pores will be used in the following calculations. Next we generate a geometry for the network and the phase, in this case air. A geometry can span over a part of the network only, so we need to specify to which pores and throats this geometry object should apply. For this example, we want it to apply to all pores and throats of the network. To do so, we can get all pore and throat indices of the network with the ``pn.pores()`` and ``pn.throats()`` calls, and pass these to the geometry object.
	
.. code-block:: python

	Ps = pn.pores()
	Ts = pn.throats()
	geo = OpenPNM.Geometry.Stick_and_Ball(network=pn,name='basic',pores=Ps,throats=Ts)
	air = OpenPNM.Phases.Air(network=pn)
	water = OpenPNM.Phases.Water(network=pn)
    
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Define the Pore-scale Physics
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

To perform most algorithms, it is necessary to define the pore scale physics that relates pores/throat geometry, phase properties, and the mechanism to be modeled.  In the case of drainage curve simulation, it is necessary to define the pressure at which the invading phase can enter each throat.  This is very commonly done by assuming the throat is a cylinder and using the so-called 'Washburn' equation.  The OpenPNM Physics module has a submodule for capillary pressure methods, including the Washburn model.  To use this model in a simulation, you first create a generic Physics object. As physics objects can span over several geometries, we need to specify to which pores and throats this physics object should apply. In this case, we want it to apply to all pores and throats. This can be done by getting all the pores and throats indices of the network with the ``pn.pores()`` and ``pn.throats()`` calls, and pass these to the physics object. 

.. code-block:: python
	
	phys = OpenPNM.Physics.Standard(network=pn,phase=water,pores=Ps,throats=Ts)

Then add the desired methods to this object using:

.. code-block:: python

    phys.add_model(propname='throat.pc',model=OpenPNM.Physics.models.capillary_pressure.washburn)


This means that the Physics object will now have a function called 'capillary_pressure', that when called will calculate throat entry pressures using the 'washburn' model.  The Washburn model requires that the Phase object (Water in this case) has the necessary physical properties of surface tension and contact angle.  

.. note::

	Both surface tension and contact angle are actually 'phase system' properties, rather than solely water properties.  It is an open problem in OpenPNM to figure out how to treat these sort of properties more rigorously.  For the present time, they must be entered a single phase properties.
	
The predefined Water object is assigned a contact angle of 110 degrees by default (water on Teflon). To change this value, it is simply a matter of replacing the 'contact_angle' function attached to the Water object with a new function that calculates the desired contact angle, as shown below:


.. code-block:: python

    water.add_model(propname='pore.contact_angle',model=OpenPNM.Phases.models.misc.constant,value=140)

In this case, the contact angle is not so much calculated as assigned a value of 140, but the result is the same.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Run a Drainage Simulation
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

At this point, the system is fully defined and ready to perform some simulations.  A typical algorithm used in pore network modeling is to use ordinary percolation to simulate drainage of wetting phase by invasion of a nonwetting phase.  An Algorithm object is be created as follows:

.. code-block:: python

	OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=pn,loglevel=20)

Before performing simulations with this algorithm it is necessary to specify the desired experimental parameters in the ``run()`` command:

.. code-block:: python
	
	Ps = pn.pores(labels=['bottom_face'])
	OP_1.run(invading_phase=water,defending_phase=air,inlets=Ps)
	
The first line in the above block finds all the pores in the network that are labeled 'bottom_face' and assigns it to 'Ps'.  This labeling step was applied during the network construction.  The list of pores which are to be considered as phase inlets along with which phases are the invading and defending phase are set to the `run()` method and the algorithm proceeds.  Upon completion one can view resulting capillary pressure curving using the following command:

.. code-block:: python

	OP_1.plot_drainage_curve()

-------------------------------------------------------------------------------
Sharing Algorithm Results Throughout the Simulation
-------------------------------------------------------------------------------

The results of the above simulation (and all simulations) are stored locally on the algorithm object.  If these results are to be used in other parts of the simulations, then they must be explicitly sent 'out'.  Keeping the results *silo-ed* in this way prevents unintentional overwriting of results by subsequent algorithms.  This allows for multiple simulations of the same type to be run with different conditions and such.  Sending the results of any simulation 'out' is done by with the `update()` command:

.. code-block:: python
	
	OP_1.update(Pc=8000)

The above command outputs data called 'occupancy' to the invading phase object. This data describes which pores and throats are filled by invading and defending phase at the specified applied capillary pressure *Pc*.  This information can be used by subsequent algorithms.  For instance it is often of interest to determine the gas phase diffusivity through a partially water filled network.  The Fickian diffusion algorithm then would use this information and set gas diffusion through water filled pores to zero and a relative effective diffusivity value could be found. 