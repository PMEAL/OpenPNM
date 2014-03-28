.. _drainage-example:

===============================================================================
Drainage Curve on a Cubic Network
===============================================================================

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Generating the Network, adding Geometry and creating Fluids
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Start by generating a basic cubic network and the other required components:

.. code-block:: python

    import OpenPNM
    pn = OpenPNM.Network.Cubic(name='test').generate(lattice_spacing=[0.0001],divisions=[10,10,10],add_boundaries=True)
    geo = OpenPNM.Geometry.Stick_and_Ball(network=pn,name='basic')
    geo.regenerate()
    air = OpenPNM.Fluids.Air(network=pn)
    air.regenerate()
	water = OpenPNM.Fluids.Water(network=pn)
	water.regenerate()

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Define the Pore-scale Physics
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

To perform most algorithms, it is necessary to define the pore scale physics that relates pores/throat geometry, fluid properties, and the mechanism to be modeled.  In the case of drainage curve simulation, it is necessary to define the pressure at which the invading fluid can enter each throat.  This is very commonly done by assuming the throat is a cylinder and using the so-called 'Washburn' equation.  The OpenPNM Physics module has a submodule for capillary pressure methods, including the Washburn model.  To use this model in a simulation, you first create a generic (empty) Physics object.  

.. code-block:: python
	
    phys = OpenPNM.Physics.GenericPhysics(network=pn,fluid=water,geometry=geo,name='phys')

Then add the desired methods to this object using:

.. code-block:: python

    phys.add_method(prop='capillary_pressure',model='washburn')
    phys.regenerate()

This means that the Physics object will now have a function called 'capillary_pressure', that when called will calculate throat entry pressures using the 'washburn' model.  The Washburn model requires that the Fluid object (Water in this case) have the necessary physical properties of surface tension and contact angle.  

.. note::

	Both surface tension and contact angle are actually 'fluid system' properties, rather than soley water properties.  It is an open problem in OpenPNM to figure out how to treat these sort of properties more rigorously.  For the present time, they must be entered a single phase properties.
	
The predefined Water object is assigned a contact angle of 110 degrees by default (water on Teflon). To change this value, it is simply a matter of replacing the 'contact_angle' function attached to the Water object with a new function that calculates the desired contact angle, as shown below:


.. code-block:: python

    water.add_method(prop='contact_angle',model='constant',value=140)
    water.regenerate()

In this case, the contact angle is not so much calculated as assigned a value of 140, but the result is the same.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Run a Drainage Simulation
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

At this point, the system is fully defined and ready to perform some simulations.  A typical algorithm used in pore network modeling is to use ordinary percolation to simulate drainage of wetting phase by invasion of a nonwetting phase.  An Algorithm object is be created as follows:

.. code-block:: python

	OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=pn, name='OP_1')

Before performing simulations with this algorithm it is necessary to setup the desired experimental parameters using ``setup()``, and then call ``run()``:

.. code-block:: python
	
	injection_sites = pn.get_pore_indices(labels='bottom')
	OP_1.setup(invading_fluid='water',defending_fluid='air',inlets=injection_sites,npts=20)
	OP_1.run()
	
The first line in the above block finds all the pores in the network that are labeled 'bottom'.  This labeling step was applied during the network construction.  The list of pores which are to be considered as fluid inlets along with which fluids are the invader and defender are set to the `run()` method and the algorithm proceeds.  Upon completion one can view resultant capillary pressure curving using `OP_1.plot_drainage_curve()`.

-------------------------------------------------------------------------------
Sharing Algorithm Results Throughout the Simulation
-------------------------------------------------------------------------------

The results of the above simulation (and all simulations) are stored locally on the algorithm object.  If these results are to be used in other parts of the simulations, then they must be explicitly sent 'out'.  Keeping the results *silo-ed* in this way prevents unintentional overwriting of results by subsequent algorithms.  This allows for multiple simulations of the same type to be run with different conditions and such.  Sending the results of any simulation 'out' is done by with the `update()` command.  Each algorithm :

.. code-block:: python
	
	OP_1.update(Pc=8000)

The above command outputs data called 'occupancy' to the invading fluid object. This data describes which pores and throats are filled by invading and defending fluid at an applied capillary pressure of 5000.  This information can be used by subsequent algorithms.  For instance it is often of interest to determine the gas phase diffusivity through a partially water filled network.  The Fickian diffusion algorithm then would use this information and set gas diffusion through water filled pores to zero and a relative effective diffusivity value could be found. 


 




















