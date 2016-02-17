.. _intermediate_usage:

###############################################################################
Tutorial 2 of 3: Digging Deeper with OpenPNM
###############################################################################

This tutorial will follow the same outline as the :ref:`getting started tutorial <getting_started>`, but will dig a little bit deeper at each step to reveal the more advanced features and usage of OpenPNM.  It is recommend that you complete that tutorial before attempting this one.

As usual, start by importing the OpenPNM package, along with the Scipy package which is always handy:

>>> import OpenPNM
>>> import scipy as sp

===============================================================================
Building a Cubic Network
===============================================================================

Let's generate a cubic network again, but with a different connectivity:

>>> pn = OpenPNM.Network.Cubic(shape=[20, 20, 10],
...                            spacing=0.0001,
...                            connectivity=8)

This **Network** has pores distributed in a cubic lattice, but connected to diagonal neighbors due to the ``connectivity`` being set to 8 (the default is 6 which is othogonal neighbors).  The various options are outlined in the *Cubic* class's documentation which can be viewed with the Object Inspector in Spyder.  OpenPNM includes several other classes for generating networks including random topology based on Delaunay tessellations (**Delaunay**).  It is also possible to import networks from external code that extracts networks from tomographic images (see [refs]).

===============================================================================
Initialize and Build Geometry Objects
===============================================================================

In this tutorial we will make a material that has different geometrical properties in different regions.  This will demonstrate the motivation behind separating the **Geometry** properties from the **Network** topology.  Let's say that the pores on the top and bottom surfaces are smaller than the internal pores.  We need to create one **Geometry** object to manage the top and bottom pores, and a second to manage the remaining internal pores:

>>> Ps1 = pn.pores(['top', 'bottom'])
>>> Ts1 = pn.find_neighbor_throats(pores=Ps1, mode='union')
>>> geom1 = OpenPNM.Geometry.GenericGeometry(network=pn, pores=Ps1, throats=Ts1,
...                                          name='surface')
>>> Ps2 = pn.pores(['top', 'bottom'], mode='not')
>>> Ts2 = pn.find_neighbor_throats(pores=Ps2, mode='intersection')
>>> geom2 = OpenPNM.Geometry.GenericGeometry(network=pn, pores=Ps2, throats=Ts2,
...                                          name='core')

The above statements result in two distinct **Geometry** objects, each applying to different regions of the domain. ``geom1`` applies to only the pores on the top and bottom surfaces (atuomatically labeled 'top' and 'bottom' during the network generation step), while ``geom2`` applies to the pores 'not' on the top and bottom surfaces.

The assignment of throats is more complicated and also illustrates the ``find_neighbor_throats`` method, which is one of the more useful topological query methods on the **Network** class.  In both of these calls, all throats connected to the given set of pores (``Ps1`` or ``Ps2``) are found; however, the ``mode`` argument alters which throats are returned.  The terms 'union' and 'intersection' are used in the *set theory* sense, such that 'union' returns *all* throats connected to the pores in the supplied list, while 'intersection' returns the throats that are *only* connected to the supplied pores.  More specifically, if pores 1 and 2 have throats [1, 2] and [2, 3] as neighbors, respectively, then the 'union' mode returns [1, 2, 3] and the 'intersection' mode returns [2].

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Add Desired Properties to Geometry
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

In the :ref:`getting started tutorial <getting_started>` we only assigned static values to **Geometry** object, which we calculated explicitly.  In this tutorial we will use the *pore-scale models* that are provided with OpenPNM.

Before applying models, however, let's assign a static random seed value between 0 and 1 to each pore on both **Geometry** objects.  We will then use these seed values in pore-scale models to generate actual pores diameters from statistical distribution functions.  To create the small pores on the surface of the domain we will adjust the parameters used in the statistical distribution.  The need to maintain two distinct sets of parameters is the driving force for defining two **Geometries**.  To start, let's put random numbers into each Geometry's 'pore.seed' property:

>>> geom1['pore.seed'] = sp.rand(geom1.Np)
>>> geom2['pore.seed'] = sp.rand(geom2.Np)

It is crucial to note that the above lines each produced an array of different length, corresponding to the number of pores assigned to each **Geometry** object.  This is accomplished by the calls to ``geom1.Np`` and ``geom2.Np``, which return the number of pores on each object.  Every Core object in OpenPNM possesses the same set of methods for mananging their data, such as counting the number of pore and throat values they represent; thus, ``pn.Np`` returns 1000 while ``geom1.Np`` and ``geom2.Np`` return 200 and 800 respectively.  The segmentation of the data between separate Geometry objects is essential to the management of pore-scale models, as will be explained next.

-------------------------------------------------------------------------------
Add Pore Size Distribution Models to Each Geometry
-------------------------------------------------------------------------------

Pore-scale models are mathematical functions that are applied to each pore (or throat) in the network to produce some local property value.  Each of the modules in OpenPNM (Network, Geometry, Phase and Physics) have a "library" of pre-written models located under "models" (i.e. *Geometry.models*).  Below this level, the models are further categorized according to what property they calculate, and there are typical 2-3 models for each.  For instance, under ``Geometry.models.pore_diameter`` you will see ``random``, ``normal`` and ``weibull`` among others.

Pore size distribution models are assigned to each Geometry object as follows:

>>> geom1.models.add(propname='pore.diameter',
...                  model=OpenPNM.Geometry.models.pore_diameter.normal,
...                  scale=0.00002, loc=0.000001,
...                  seeds = 'pore.seed')
>>> geom2.models.add(propname='pore.diameter',
...                  model=OpenPNM.Geometry.models.pore_diameter.weibull,
...                  shape=1.2, scale=0.00005, loc=0.000001,
...                  seeds = 'pore.seed')

Pore-scale models tend to be the most complex (and confusing) aspects of OpenPNM, so it's worth dwelling on the above two commands:

(1) Both ``geom1`` and ``geom2`` have a ``models`` attribute where the parameters specified in the ``add`` command are stored for future use if/when needed.  The ``models`` attribute actually contains a **ModelsDict** object which is a customized dictionary for storing and managing this type of information.  Details of the **ModelsDict** class are outlined elsewhere [???].

(2) The ``propname`` argument specifies which property the model calculates.  This means that the numerical results of the model calculation will be saved in their respective **Geometry** objects as ``geom1['pore.diameter']`` and ``geom2['pore.diameter']``.

(3) Each model stores it's result under the same ``propname`` but these values do not conflict since each **Geometry** object presides over a unique set of pores and throats.

(4) The ``model`` argument contains a *handle* to the desired function, which is extracted from the *models* library of the relevant *Module* (**Geometry** in this case).  Each **Geometry** object has been assigned a different statistical model, *normal* and *weibull*.  This ability to apply different models to different regions of the domain is reason multiple **Geometry** objects are permitted.  The added complexity is well worth the added flexibility.

(5) The remaining arguments are those required by the chosen *model*.  In the above cases, these are the parameters that define the statistical distribution.  Note that the mean pore size for ``geom1`` will be 20 um (set by ``scale``) while for ``geom2`` it will be 50 um, thus creating the smaller surface pores as intended.  The pore-scale models are well documented regarding what arguments are required and their meaning; as usual these can be viewed with Object Inspector in Spyder.

-------------------------------------------------------------------------------
Add Additional Pore-Scale Models to Each Geometry
-------------------------------------------------------------------------------

In addition to pore diameter, there are several other geometrical properties needed to perform a permeability simulation.  Let's start with throat diameter:

>>> geom1.models.add(propname='throat.diameter',
...                  model=OpenPNM.Geometry.models.throat_misc.neighbor,
...                  pore_prop='pore.diameter',
...                  mode='min')
>>> geom2.models.add(propname='throat.diameter',
...                  model=OpenPNM.Geometry.models.throat_misc.neighbor,
...                  pore_prop='pore.diameter',
...                  mode='min')

Instead of using statistical distribution functions, the above lines use the ``neighbor`` model which assigns each throat the value of the specified 'pore_prop' from it's neighboring pores.  In this case, each throat is assigned the minimum pore diameter of it's two neighboring pores.  Other options for ``mode`` include 'max' and 'mean'.

We'll also need throat length as well as the cross-sectional area of pores and throats, for calculating the hydraulic conductance model later.

>>> geom1.models.add(propname='throat.length',
...                  model=OpenPNM.Geometry.models.throat_misc.straight)
>>> geom2.models.add(propname='throat.length',
...                  model=OpenPNM.Geometry.models.throat_misc.straight)
>>> geom1.models.add(propname='throat.area',
...                  model=OpenPNM.Geometry.models.throat_misc.cylinder)
>>> geom2.models.add(propname='throat.area',
...                  model=OpenPNM.Geometry.models.throat_misc.cylinder)
>>> geom1.models.add(propname='pore.area',
...                  model=OpenPNM.Geometry.models.throat_misc.sphere)
>>> geom2.models.add(propname='pore.area',
...                  model=OpenPNM.Geometry.models.throat_misc.sphere)

At this point you might ask "'"*why can't we just calculate pore and throat cross-sectional areas manually and assign them as we did in* :ref:`Tutorial #1 <getting_started>`"?  The answer is 'you can, but you shouldn't'.  The reason is that pore-scale models can be 'recalculated' or 'regenerated', so changes in one property can be easily reflected in all dependent properies.  For instance, if you wish to perform a simulation on a new realization of the network, you only need to alter the random seed values assigned to ``geom1`` and ``geom2`` the 'regenerate' all the models as follows:

>>> geom1['pore.seed'] = sp.rand(geom1.Np)
>>> geom2['pore.seed'] = sp.rand(geom2.Np)
>>> geom1.models.regenerate()
>>> geom2.models.regenerate()

The first two lines assign new random numbers to each pore, and the final two lines cause all of the pore-scale models to be recalculated, using the same parameters specified above.  This means that all pore diameters change (but still following the same statistical distribution), thus so will the throat diameters which were taken as the minimum of the two neighboring pores, and so on.  Note that during the regeneration process all models are called in the order they were originally added.

===============================================================================
Initialize and Build Phase Objects
===============================================================================

**Phase** objects are defined in similar manner to the **Geometry** objects outlined above.  For this tutorial, we will create two generic **Phase** objects for 'air' and 'water', and assign pore-scale models for calculating viscosity which is the only phyical property needed for permeability calculations.

>>> air = OpenPNM.Phases.GenericPhase(network=pn)
>>> water = OpenPNM.Phases.GenericPhase(network=pn)

A variety of pore-scale models are available for calculating **Phase** properties.  These models are generally taken from correlations in the literature:

>>> air.models.add(propname='pore.viscosity',
...                model=OpenPNM.Phases.models.viscosity.chung)
>>> water.models.add(propname='pore.viscosity',
...                  model=OpenPNM.Phases.models.viscosity.water)

Note that all **Phase** objects are automatically assigned standard temperature and pressure conditions when created.  This can be adjusted:

>>> air['pore.temperature'] = 333  # K
>>> water['pore.temperature'] = 353  # K

Since viscosity is highly dependent on temperature, it is necessary to 'regeneate' the viscosity models:

>>> air.models.regenerate()
>>> water.models.regenerate()

===============================================================================
Initialize and Build Physics Objects
===============================================================================
