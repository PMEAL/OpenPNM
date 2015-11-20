.. _IP-example:

===============================================================================
Predicting Fluid Configurations using Invasion Percolation
===============================================================================

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Generating the Network, adding Geometry and creating Phases
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Start by generating a basic cubic network and the other required components.  We'll use a 2D network to illustrate the algorithm since it's easier for visualization.

.. code-block:: python

	import OpenPNM
	pn = OpenPNM.Network.Cubic(shape=[20, 20, 1], spacing=0.0001)

Next we need to create a Geometry object to manage the pore and throat size information, and a phase object to manage the thermophysical properties of the invading fluid.

.. code-block:: python

	geom = OpenPNM.Geometry.Toray090(network=pn,
	                                 pores=pn.Ps,
																	 throats=pn.Ts)
	water = OpenPNM.Phases.Water(network=pn)

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Define the Pore-scale Physics
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

To perform most algorithms, it is necessary to define the pore scale physics that relates pore/throat geometry, phase properties, and the mechanism to be modeled.  In the case of a drainage curve simulation, it is necessary to define the pressure at which the invading phase can enter each throat.  This is commonly done by assuming the throat is a cylinder and using the so-called 'Washburn' equation.  The OpenPNM Physics module has a submodule for capillary pressure methods, including the Washburn model.  To use this model in a simulation, you first create a generic Physics object, then add the desired methods to this object as follows:

.. code-block:: python
	phys = OpenPNM.Physics.GenericPhysics(network=pn,
                                        phase=water,
																				geometry=geom)
	phys.add_model(propname='throat.capillary_pressure',
	               model=OpenPNM.Physics.models.capillary_pressure.washburn)

This means that the Physics object will now have a model called 'throat.capillary_pressure', that when called will calculate throat entry pressures using the 'washburn' model.  For more on the use of models see :ref:`Pore Scale Models<models>`.  The Washburn model requires that the Phase object (Water in this case) has the necessary physical properties of surface tension and contact angle.

The predefined Water object is assigned a contact angle of 110 degrees by default (water on Teflon). This can be confirmed by printing the first element of the 'pore.contact_angle' array:

>>> print(water['pore.contact_angle'])[0]
110.0

To change this value, the 'pore.contact_angle' property of the Water object can be set to a new constant:

>>> water['pore.contact_angle'] = 140.0

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Run an Invasion Simulation
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

At this point, the system is fully defined and ready to perform some simulations.  An IP algorithm object is instantiated as follows:

>>> IP = OpenPNM.Algorithms.InvasionPercolation(network=pn, phase=water)

Before running the algorithm it is necessary to specify the inlet sites from where the invading fluid enters the network:

>>> IP.set_inlets(pores=pn.pores('left'))

The final step is to invaded the network.  This is accomplished with the ``run`` method of the IP object.

>> IP.run()

This method produces arrays called 'pore.invaded' and 'throat.invaded' on the IP object that contain the invasion sequence of each pore and throat, respectively.  To obtain the invading fluid configuration at some intermediate invasion state (for instance the first 200 invasions) it is simply a matter of applying a boolean operated to the 'pore.invaded' and 'throat.invaded' arrays such as:

>>> Pinv = IP['pore.invaded'] < 200
>>> Tinv = IP['throat.invaded'] < 200

These partially invaded configurations could then be used to perform various simulations such as relative permeability.

The ``run`` command takes an option argument of ``n_steps``, which if given performs a partial invasion of the network.  This approach is requried if you wish to perform more complex invasions such as 100 steps from the 'left', then 100 from the 'right'.  In this case the invasion pattern will not be the same as if the invasion had proceeded entirely from the 'left'.

.. note:: Advanced IP

    A more sophisticated version of the IP algorithm is currently under development that accepts the volume of injection rather than just the number of steps.  It also allows for multiple clusters be invaded at the same time, and for each defending cluster to be invaded at a different rate.  The basis IP algorithm described above will remain in the code since it is handy to have simple algorithm in many cases.
