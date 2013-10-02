===============================================================================
Pore Scale Physicis
===============================================================================

-------------------------------------------------------------------------------
General Usage
-------------------------------------------------------------------------------
In OpenPNM each physical process is given its own class(? - file), and within each class any number of models/equations can be defined by simply adding the appropriate method/equations.  This setup was designed to facilitate easy editability of the physics models which can be highly varied.  

-------------------------------------------------------------------------------
Capillary Pressure
-------------------------------------------------------------------------------

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Washburn
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The Washburn equation, often refered to as the Young-Laplace equation, is the most commonly used means to relate throat size to capillary entry pressure.

.. math::

  P_c = {\frac{-4 \sigma cos(\theta)}{D}}

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Purcell Toroid
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The Purcell toroid is a more sophisticated means of relating throat size to capillary entry pressure that accounts for the converging-diverging nature of throats.  This model can lead to significantly increased entry pressures, especially for neutrally wetting systems.  


-------------------------------------------------------------------------------
Flow
-------------------------------------------------------------------------------

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Hagen-Poiseuille
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


-------------------------------------------------------------------------------
Diffusion
-------------------------------------------------------------------------------

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Knudsen
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bulk
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


-------------------------------------------------------------------------------
Writing Custom Physics
-------------------------------------------------------------------------------
To add a model to an existing physics class, simple add a new def with the desired equation.

To add a new physics class, create a new file, and notify OpenPNM to open it by add it to the Physics/__init.py__.
