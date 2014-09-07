.. _io:

===============================================================================
Saving and Loading Objects and Entire Simulations
===============================================================================
OpenPNM has the ability to save and and load both individual objects and full simulations.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Individual Objects
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The OpenPNM.Base class contains a ``save`` method, which means that it is inherited by all the main OpenPNM objects.  This method is quite simplistic, and only saves the objects data (i.e. *properties* and *labels*).  It does not save any models on the object or any information about the simulation hierarchy.  This method is really just to save the pure data arising from a simulation, or time consuming geometric calculation.  The ``load`` method that is also part of the Base class is able to read in files created by ``save`` and populate an empty Generic object with the data.  Loading the data onto a Generic object allows the data to be viewed and manipulated using the relevant OpenPNM methods.  The following code block demonstrates this:

.. code-block:: python

    import OpenPNM
	pn = OpenPNM.Network.Cubic(shape=[3,3,3])
	pn.save('test_pn')
	
	#Now create and empty generic Network
	gn = OpenPNM.Network.GenericNetwork()
	#Use the load method to retrieve the saved data and place into the empty Network object
	gn.load('test_pn')
	
In addition to loading the data into the dictionary of the generic object, the ``load`` method also changes the name of the 'loading' object.  This means that the new object name will match any label names that may have been created by the old object.  This would be essential for trying to piece together a simulation from saved objects.  Of course, if saving a simulation is the aim, then the IO module in OpenPNM.Utilities is a better option.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Full Simulations
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The IO module located under OpenPNM.Utilities contains several classes for saving more then just individual objects.  The PNM class is designed specifically for saving a simulation in it's entirety, so that it can be reloaded and used for further simulations.  

.. code-block:: python

    import OpenPNM
	pn = OpenPNM.Network.Cubic(shape=[3,3,3])
	geo = OpenPNM.Geometry.Stick_and_Ball(network=pn,pores=pn.pores(),throats=pn.throats(),name='geo_1')
	air = OpenPNM.Phases.Air(network=pn)
	phys = OpenPNM.Physics.Standard(network=pn,phase=air,pores=pn.pores(),throats=pn.throats())
	
	import OpenPNM.Utilities.IO as io
	io.PNM.save(pn,'test_pn')
	
The ``PNM.save`` creates a specialized '.pnm' file that contains all the necessary information to recreate the simulation.  It can be reloaded with:

.. code-block:: python

    import OpenPNM
    import OpenPNM.Utilities.IO as io
    pn = io.PNM.load('test_pn')

This procedure returns a network object only, but the network retains a link to all the objects with which is was associated before being saved.  These links can be accessed using the ``geometries`` , ``phases`` and ``physics`` methods.  For instance, to obtain a handle to the 'geo' object:

>>> geo = pn.geometries('geo_1')


.. warning:: 
    
	There is currently an important limitation on the PNM save/load features: it does not retain the class type of the saved object.  This is acceptable for the Geometry, Phase and Physics objects, but most Network objects have additional methods added (such as ``asarray`` and ``fromarray``).  These methods would not be available to the loaded object.  

The IO module also includes the ability to output to VTK and Matlab MAT files.  


