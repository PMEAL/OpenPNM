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
|                | float      | of pores to be imported          |
+----------------+------------+----------------------------------+
| pdiameter      | <Npx1>     | pore diamters, in meters         |
|                | float      |                                  |
+----------------+------------+----------------------------------+
| pvolume        | <Npx1>     | pore volumes, in cubic meters    |
|                | float      |                                  |
+----------------+------------+----------------------------------+
| pnumbering     | <Npx1>     | = 0:1:Np-1                       |
|                | int        |                                  |
+----------------+------------+----------------------------------+
| tconnections   | <Ntx2>     | pore numbers of the two pores    |
|                | int        | that each throat connects        |
+----------------+------------+----------------------------------+
| tdiameter      | <Ntx1>     | throat diameters, in meters      |
|                | float      |                                  |
+----------------+------------+----------------------------------+
| tnumbering     | <Ntx1>     | = 0:1:Nt-1                       |
|                | int        |                                  |
+----------------+------------+----------------------------------+

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Importing with MAT
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Once you have correctly formatted a *.mat file, it can be loaded with the following commands.

.. code-block:: python
    
    import OpenPNM
    fname = 'examples/yourfile' # or 'examples/yourfile.mat'
    pn = OpenPNM.Network.MatFile(name='mat_net',filename=fname)
    geom = pn.gemetries('internal')

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Additional Pore and Throat Properties
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Additional properties and labels can be added to the network as long as they are arrays of the correct length. For example, if you saved the *.mat file with the variable `pshape` that was an Npx1 float array, you could include it in your import as follows:

>>> pn = OpenPNM.Network.MatFile(name='mat_net',filename=fname,xtra_pore_data='shape')
>>> # this will fail unless you actually build a *.mat file with a variable named "pshape"

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Adding Surfaces and Boundaries to Network with ptype and ttype
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
There is no add_boundaries() command for this Network class. But, you can add `ptype` and `ttype` variables to the *.mat file in order to store "surface" and "boundary" information. In OpenPNM, "boundary" pores are zero-volume pores just outside the network domain, that typically gain boundary conditions in simulations. Conventionally, there is one "boundary" pore placed for every "surface" pore in the system, where "surface" pores are pores on the very edge of the domain.

Currently, this class has a built-in, but rather inflexible method that assumes the users knows what they are doing. Follow these instructions carefully:

In order for the nework to import boundaries, the pore and throat data `ptype` and `ttype` must be present. 

+----------------+------------+----------------------------------+
| Variable Name  | Value      | Description                      |
+================+============+==================================+
| ptype          | <Npx1>     | (optional) designates surfaces   |
|                | int        | of pores in network.             |
+----------------+------------+----------------------------------+
| ttype          | <Ntx1>     | (optional) designates surfaces   |
|                | int        | of throats in network.           |
+----------------+------------+----------------------------------+

The `type` variables are integers between 0 and 6. All internal pores, inlcuding "surface" pores, should have the value 0. The rest should be labelled as follows: 1-top, 2-left, 3-front, 4-back, 5-right, 6-bottom.

Importing a network with boundaries can be done as follows:

>>> pn = OpenPNM.Network.MatFile(name='mat_net',filename=fname,xtra_pore_data='type',xtra_throat_data='type')

Because the `type` variables are loaded, the importer automatically adds labels for boundaries and surfaces.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Example File
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

One example *.mat file has been loaded with the OpenPNM installation. Look for /OpenPNM/Network/test_pn.mat


