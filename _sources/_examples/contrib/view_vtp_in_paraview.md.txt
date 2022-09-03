# Using Paraview to Visualize Output

This tutorial will show step by step how to plot the networks into Paraview. This post-processing software is open-source and can be downloaded from the [Paraview Website](http://www.paraview.org/).

## Open the VTK file created by OpenPNM

Once the OpenPNM script (run_script.py) have been run, a output file named net.vtp is created. This file can be opened in Paraview and contains all the parameters computed by OpenPNM. A summary of them is seen under the Properties bar. Once the file is open, click on Apply, this action should be done after each modification.

1. Open the net.vtp file

![](http://i.imgur.com/gmPVmRM.png)

2. Click and apply, the following network should be seen:

![](http://i.imgur.com/wFinDmX.png)


## Plot the pores data

To visualize the pore data, we need to add some glyphs. First click on the Glyph button in the tool bar. Then, you can plot the pore data as sphere where their size and/or their colour is related to some variables. In the images below we set up spheres pore, where the diameter is link to the pore diameters computed by OpenPNM and the colour to the mole fraction. Click on Apply button visualize this settings.

Pore data plotted with spheres radius proportional to pore diameter and sphere colour linked to the pore mole fraction:

![](http://i.imgur.com/bnFuH3r.png)

## Plot the throat data

To visualize the throat data (like their diameters or the molar flow rates) we need to set up the following filters in Paraview. First, use the Shrink filter and set it up to 1. Then, the cell data needs to be transposed to the data point with CellDatatoPointdata filter. Then extract the surface with the filter ExtractSurface. Finally, the data may be plotted as tube by using the Tube filter. As previously for the pore data, either the Tube radius or colour can be linked to throat data.

Throat data plotted with tubes radius proportional to throat diameter and tubes colour linked to the throat mole rate:

![](http://i.imgur.com/SX5YeVj.png)

## Plot different water saturations

To visualize the water invasion sequence we need to use the threshold filter. It has to be added after the Glyph filter. In the threshold properties tool bar, set up the water.pore.inv_seq as the threshing variable and a solid colour in Colouring tool bar. The same operation may be done for the throat water invasion sequence.

Pores water invasion sequence, threshold 7:

![](http://i.imgur.com/lfowwsV.png)

Pores + Throats water invasion sequence, threshold 7:

![](http://i.imgur.com/dmmPNCW.png)
