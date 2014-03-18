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

	OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=pn,name='OP_1')

To perform simulations using this algorithm simply call the `run()` command with the desired parameters:

.. code-block:: python
	
	injection_sites = pn.get_pore_indices(labels='bottom')
	OP_1.run(invading_fluid='water',defending_fluid='air',inlets=injection_sites,npts=20)
	
The first line in the above block finds all the pores in the network that are labeled 'bottom'.  This labeling step was applied during the network construction.  The list of pores which are to be considered as fluid inlets along with which fluids are the invader and defender are set to the `run()` method and the algorithm proceeds.  Upon completion one can view resultant capillary pressure curving using `OP_1.plot_drainage_curve()`.

-------------------------------------------------------------------------------
Sharing Algorithm Results Throughout the Simulation
-------------------------------------------------------------------------------

The results of the above simulation (and all simulations) are stored locally on the algorithm object.  If these results are to be used in other parts of the simulations, then they must be explicitly sent 'out'.  Keeping the results *silo-ed* in this way prevents unintentional overwriting of results by subsequent algorithms.  This allows for multiple simulations of the same type to be run with different conditions and such.  Sending the results of any simulation 'out' is done by with the `update()` command.  Each algorithm :

.. code-block:: python
	
	OP_1.update(Pc=8000)

The above command outputs data called 'occupancy' to the invading fluid object. This data describes which pores and throats are filled by invading and defending fluid at an applied capillary pressure of 5000.  This information can be used by subsequent algorithms.  For instance it is often of interest to determine the gas phase diffusivity through a partially water filled network.  The Fickian diffusion algorithm then would use this information and set gas diffusion through water filled pores to zero and a relative effective diffusivity value could be found. 


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Run a Diffusion Simulation a Partially Saturated Network
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Calculating the gas phase diffusivity through a water invading porous medium is one of the main applications of pore networks.  The Fickian diffusion algorithm supplied with OpenPNM is setup and called in much the same way as the ordinary percolation algorithm described above.

-------------------------------------------------------------------------------
Prepare the Algorithm
-------------------------------------------------------------------------------

Firstly, an Algorithm object must be instantiated:

.. code-block:: python

	Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(network=pn,name='Fickian_alg')

Each algorithm performs drastically different functions and calculations so each naturally expect quite different arguments.  The Fickian algorithm needs to know what boundary conditions are prevailing.  These can include Dirchelet, various types of Neumann, reaction rates, and so on.  The lines below outline how to setup Dirchelet conditions on two opposing faces.  Note that the process involves first finding the pore finding the indices of pores laying on the 'top' or 'bottom' face of the domain, then applying a the 'Dirichlet' label, and finally applying the boundary value to those locations. 

.. code-block:: python

	BC1 = pn.get_pore_indices(labels=['top'],mode='intersection')
	Fickian_alg.set_pore_info(label='Dirichlet', locations=BC1)
	Fickian_alg.set_pore_data(prop='BCval', data=0.6, locations=BC1)
	BC2 = pn.get_pore_indices(labels=['bottom'],mode='intersection')
	Fickian_alg.set_pore_info(label='Dirichlet', locations=BC2)
	Fickian_alg.set_pore_data(prop='BCval', data=0.2, locations=BC2)
	
Note that this simulation will run on a Network that has been invaded upto 5000 Pa with water due to the OP_1.update(Pc=5000) command used above.  It is a simple matter to change the network saturation be calling this command with a different applied pressure.  

There are many features, details and nuances of this package that have been glossed over in this quickstart guide.  The complete documentation describes the OpenPNM framework in detail.  Happy coding.  




















