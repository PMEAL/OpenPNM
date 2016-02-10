.. _intermediate_usage:

###############################################################################
Tutorial 2 of 3: Digging Deeper with OpenPNM
###############################################################################

This tutorial will follow the same outline as the :ref:`getting started tutorial <getting_started>`, but will dig a little bit deeper at each step to reveal the more advanced features and usages of OpenPNM.  Be sure you've done and understood that tutorial before attempting this one.

As usual, start by importing the OpenPNM package and the Scipy package which is always handy:

>>> import OpenPNM
>>> import scipy as sp

===============================================================================
Building a Cubic Network
===============================================================================

Let's generate a cubic network but with a different connectivity:

>>> pn = OpenPNM.Network.Cubic(shape=[20, 20, 10],
...                            spacing=0.0001,
...                            connectivity=8)

This **Network** has pores distributed in a cubic lattice, but connected to diagonal neighbors due to the ``connectivity`` being set to 8 (the default is 6).  The various options are outlined in the *Cubic* class's documentation which can be viewed with the object inspector in Spyder.  OpenPNM includes several other classes for generating networks including random topology based on Delaunay tessellations (**Delaunay**).  It is also possible to import networks from external code that extracts networks from tomographic images (see [refs]).

===============================================================================
Initialize and Build a Geometry Object
===============================================================================

In this tutorial we will make a material that has different geometrical properties in different regions.  This will demonstrate the motivation behind separating the **Geometry** properties from the **Network** topology.  Let's say that the pores on the top and bottom surfaces are smaller than the internal pores.  We need to create one **Geometry** object to manage the top and bottom pores, and a second to manage the remaining internal pores:

>>> Ps = pn.pores(['top', 'bottom'])
>>> geom1 = OpenPNM.Geometry.GenericGeometry(network=pn, pores=Ps,
...                                          name='surface')
>>> Ps = pn.pores(['top', 'bottom'], mode='not')
>>> Ts = pn.throats()
>>> geom2 = OpenPNM.Geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts,
...                                          name='core')

The above statements result in two distinct **Geometry** objects each applying to different regions of the full network domain.  As we shall see, this simplifies the data management in some important ways. ``geom1`` applies to only the pores on the top and bottom surfaces, while ``geom2`` applies to the pores 'not' on the top and bottom surfaces, as well as all the throats.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Add Desired Properties to Geometry
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

In the :ref:`getting started tutorial <getting_started>` we only assigned static values to Geometry object, which we calculated explicitly.  In this tutorial we will use the pre-written *pore-scale models* that are provide with OpenPNM.

.. note:: Pore Scale Models Explained - Level 1

	Pore-scale models are mathematical functions that are applied to each pore (or throat) in the network to produce some local property or value.  Each of the modules in OpenPNM (Network, Geometry, Phase and Physics) have a "library" of pre-written models located under "models" (i.e. Geometry.models).  Below this level, the models are further categorized according to what property they calculate, and there are typical 2-3 models for each.  For instance, under ``Geometry.models.pore_seed`` you will see ``random`` and ``spatially_correlated``.  Each of these produces a random number for each pore, but the ``spatially_correlated`` model places numbers of similar size in neighboring pores (using a convolution filter) while the ``random`` model just places numbers (just using the ``scipy.rand`` function).

For both **Geometry** objects, we will assign each pore a static random seed value between 0 and 1, and then will use these seed values in statistical distribution functions to generate actual pores diameters.  To create the small surface pores, we will  adjust the parameters used in the statistical distribution.  The need to maintain two distinct sets of parameters is the driving force for defining two **Geometries**.  To start, let's put random numbers into each Geometry's 'pore.seed' property:

>>> geom1['pore.seed'] = sp.rand(geom1.Np)
>>> geom2['pore.seed'] = sp.rand(geom2.Np)

It is crucial to note that the above lines each produced an array of different length, corresponding to the number of pores assigned to each **Geometry** object.  This is accomplished by the calls to ``geom1.Np`` and ``geom2.Np``, which return the number of pores on each.  Every Core object in OpenPNM possesses the same set of methods for mananging their data, such as counting the number of pore and throat values they represent; thus, ``pn.Np`` returns 1000 while ``geom1.Np`` and ``geom2.Np`` return 200 and 800 respectively.  The segmentation of the data between separate Geometry objects is essential management of pore-scale models, as will be explained next:
