.. _matload:

===============================================================================
Loading Networks Saved in MATLAB
===============================================================================
OpenPNM has the ability to load networks generated in MATLAB. Saved as a specially formated *.mat file.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MAT File Format
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
If you have created a pore network in MATLAB and you would like to import it into OpenPNM. Create the following variables:

+----------------+------------+----------------------------------+
| Variable Name  | Value      | Description                      |
+================+============+==================================+
| pcoords        | <Npx3>     | physical coordinates, in meters, |
|                | double     | of pores to be imported          |
+----------------+------------+----------------------------------+
| pdiameter      | <Npx1>     | pore diamters, in meters         |
|                | double     |                                  |
+----------------+------------+----------------------------------+
| pvolume        | <Npx1>     | pore volumes, in cubic meters    |
|                | double     |                                  |
+----------------+------------+----------------------------------+
| pnumbering     | <Npx1>     | = 0:1:Np-1                       |
|                | int32      |                                  |
+----------------+------------+----------------------------------+
| ptype          | <Npx1>     | (optional) designates surfaces   |
|                | int32      | of pores in network.             |
|                |            | (more details below)             |
+----------------+------------+----------------------------------+
| tconnections   | <Ntx2>     | pore numbers of the two pores    |
|                | int32      | that each throat connects        |
+----------------+------------+----------------------------------+
| tdiameter      | <Ntx1>     | throat diameters, in meters      |
|                | double     |                                  |
+----------------+------------+----------------------------------+
| tnumbering     | <Ntx1>     | = 0:1:Nt-1                       |
|                | int32      |                                  |
+----------------+------------+----------------------------------+
| ttype          | <Ntx1>     | (optional) designates surfaces   |
|                | int32      | of throats in network.           |
|                |            | (more details below)             |
+----------------+------------+----------------------------------+

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Importing with MAT
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Once you have correctly formatted a *.mat file, it can be loaded with the following commands.

.. code-block:: python
    
    fname = 'examples/yourfile' # or 'examples/yourfile.mat'
    pn = OpenPNM.Network.MatFile(name='mat_net',filename=fname)
    geom = pn.gemetries('internal')

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Adding Surfaces and Boundaries to Network with ptype and ttype
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
There is no add_boundaries() command for this Network class. But, you can add `ptype` and `ttype` variables to the *.mat file in order to store "surface" and "boundary" information.

Currently, this is a rather inflexible method that assumes the users knows what they are doing. 


|                |            | 0-non-surface, 1-top, 2-left, 3-front   |
|                |            | 4-back, 5-right, 6-bottom        |

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


