.. _visualize:

================================================================================
Using Paraview to Visualize Output
================================================================================

This tutorial will show step by step how to plot the networks into Paraview. This post-processing software is open-source and can be downloaded from the `Paraview Website <http://www.paraview.org/>`_.


--------------------------------------------------------------------------------
Export a VTK File from OpenPNM
--------------------------------------------------------------------------------

Assume the following project has been created and you wish to visualize it:

.. code-block:: python

    >>> import scipy as sp
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
    >>> water = op.phases.Water(network=pn)
    >>> phys = op.physics.Standard(network=pn, phase=water, geometry=geo)

Now we can export it into a variety of formats using the ``io`` module, but the
VTK format is easiest for Paraview.

.. note::

    The XDMF format is actually must better for using in Paraview if the network
    is of appreciable size since it's performance is far better, but it is
    slightly more complicated since it splits the data into 2 files.  The
    information in this tutorial applies equally well to both formats.

Now that the project is setup we can export it:

.. code-block:: python

    >>> op.io.VTK.save(network=pn, phases=water, filename='test_file')

This will create a file in your current working directory called 'test_file.vtp'.  There are several things to note about this line.  Firstly,  the 'vtp' file extension is used rather than 'vtk'. The last letter indicates the type of data in the file, where 'p' indicates PolyData.  The VTK accronym stands for Visualization Tool Kit, which is a general reference to all similar files.  Secondly, current working directory can be found from the interactive python terminal by importing the ``os`` package, and typing ``os.get_cwd()``.  Finally, sending Phase data to the VTK file is optional, and you must specify which phase data to write.

--------------------------------------------------------------------------------
Open the VTK file created by OpenPNM
--------------------------------------------------------------------------------

Once the OpenPNM script (run_script.py) have been run, a output file named net.vtp is created. This file can be opened in Paraview and contains all the parameters computed by OpenPNM. A summary of them is seen under the Properties bar. Once the file is open, click on Apply, this action should be done after each modification.

1. Open the 'test_file.vtp' file using the normal GUI functions in Paraview:

.. image:: http://i.imgur.com/gmPVmRM.png

2. Click and apply, the network should be seen in the view port looking somthing like this:

.. image:: http://i.imgur.com/wFinDmX.png

--------------------------------------------------------------------------------
Plot the Dore Data
--------------------------------------------------------------------------------

To visualize the pore data, we need to add some glyphs. First click on the Glyph button in the tool bar. Then, you can plot the pore data as spheres, where their size and/or their color can be mapped to some variables. In the images below spherical glyphs are assigned to the pores, where the diameter is linked to the pore diameters and the color to the concentration. Clicking on the Apply button renders these settings.

.. image:: http://i.imgur.com/bnFuH3r.png

--------------------------------------------------------------------------------
Plot the Throat Data
--------------------------------------------------------------------------------

.. note::

  Plotting throats as nicely oriented cylinders with controlled diameter and color remains an annoyance, although we are working on getting the necessary information into the VTK file for throats be added as simply as adding cylindrical glyphs.

To visualize the throat data (like diameters or molar flow rates) we need to set up the following filters in Paraview. First, use the Shrink filter and set it up to 1. Then, the cell data needs to be transposed to the data point with CellDatatoPointdata filter. Then extract the surface with the filter ExtractSurface. Finally, the data may be plotted as tube by using the Tube filter. As previously for the pore data, either the Tube radius or color can be linked to throat data.

Throat data plotted with tubes radius proportional to throat diameter and tubes color linked to the throat mole rate:

.. image:: http://i.imgur.com/SX5YeVj.png
