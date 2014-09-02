.. _physics:

===============================================================================
Physics
===============================================================================
Physics objects are where geometric data and fluid property data are combined to compute the pore-scale physical behavior in the simulation.  For instance, the capillary entry pressure for a throat is a function of size (from Geometry) and the surface tension of the liquid (from Phase), but there are many ways to compute the actual entry pressure, including the Washburn equation for a cylinder, or the Purcell equation for a toroid.  Specifying unique pore-scale Physics models is what sets pore network simulations apart from each other.  The Physics object manages these pore-scale properties and models.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Basic Usage
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Customizing Physics
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For description of how to create customized subclasses, add properties to the model library, and add new models see :ref:`Customizing OpenPNM<customizing>`