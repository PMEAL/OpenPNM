.. _overview:

###############################################################################
Overview of the OpenPNM Framework
###############################################################################

===============================================================================
Main Modules
===============================================================================

The OpenPNM framework is build upon 5 main objects.  

1 `Network`_: The Network object is the main 'controller' object for framework.  It contains all the topological information about the Network, as well as numerous methods for querying and manipulating the topology.  Each simulation will have only 1 Network object.

3 `Geometry`_: Geometry objects control the pore-scale geometrical properties of the network such as pore size.  A simulation may have multiple Geometry objects depending on the problem being modeled.  For instance, a stratified material may have a separate Geometry object for each layer if the pore and throat sizes differ between them.  

4 `Phases`_: Phase objects contain information about the thermo-physical properties of the liquids, gases, solids that exist in the pore space.  For instance, a Phase object for water would possess its temperature, as well as models for calculating its viscosity as a function of temperature (and any other relevant properties).

5 `Physics`_: Physics objects contain methods for calculating pore scale physics properties which combine fluid and geometry values

6 `Algorithms`_: This module is the home of the actual algorithms that use the network properties defined by the above modules


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Base
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This module contains the main abstract ``Base`` class from which all OpenPNM objects derive.  The main functionality offered by this class is the logger which outputs debugging info and error message to the console.  The ``Base`` class also contains methods for finding other objects by name, deleting objects from the simulation, and other high level operations.  

The ``Base`` module also contains the ``Core`` class, which possesses the methods for accessing and manipulating data within the framework.  All objects in OpenPNM are subclasses of the Python native ``dict`` or ``dictionary`` class.  This class is similar to a ``struct`` in other languages where data can be stored in a named, nested tree structure.  OpenPNM adds numerous methods to the default behavior of ``dict`` for doing things like labeling pores, looking them up by label, listing defined properties, etc.  ``Core`` also overrides the normal behavior of the ``dict`` 'set' and 'get' methods to help 'protect' the integrity of the data.  For instance, dictionary keys 'must' start with either 'pore' or 'throat', which indicates the type of information stored.    

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Network
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This module is the heart of OpenPNM.  It contains a ``GenericNetwork`` class which possesses suite of network query methods, based on adjacency and incidence matrices, and other graph theory concepts.  The methods in ``GenericNetwork`` are fully agnostic to the type of network due the generalized way that OpenPNM stores data.  This is explained in more detail in the :ref:`Network <network>` documentation.

This module also contains numerous subclasses of the ``GenericNetwork``, which possess the code for actually generating specific network topologies (e.g. cubic, random, etc).  All subclasses derive from GenericNetwork so have those methods, as well as any additional methods relevant to the specific topology.  For instance, the Cubic class has the ability to output a 3D array containing specified data  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Geometry
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The ``Geometry`` class contains methods for calculating pore scale properties such as pore diameter and throat length.  It's name should not be confused with ``Network`` which only controls the 'topology'.

``Geometry` objects are *built* by the user to contain the specific models to be used when calculating the pore and throat geometry properties.  For instance, a ``Geometry`` object can calculate 'pore volume' assuming the pore is a sphere, cuboid or any other shape.  The user can use the methods supplied with OpenPNM, or add their own.  

There can be multiple ``Geometry`` objects defined for different locations in a ``Network`` simultaneously.  This was intended to allow for multi-layer media (such fuel cell gas diffusion layers with microporous layers on one side), but is also quite useful when applying boundary pores which usually need to have special pore geometry such as 0 volume to produce consistent results.

The ``Geometry`` module contains a ``GenericGeometry`` class.  **Geometry** objects are essentially empty when initialized, and the user then adds the desired models to the object (a process called *composition*).  For example, it is typical to assign a random seed to each pore which is subsequently used in the calculation of pore size from some sort of statistical distribution.  OpenPNM comes with several models that can be used to generate random seeds.  These are stored in a submodule called **pore_seeds**, where the options include ``perlin_noise`` and ``random_field``.  It is envisioned that each user may add their own methods for generating seeds (or any other calculation) to this submodule (or any other), so the process of adding methods to the object was designed to be as flexible and customizable as possible.  For more information on this process see the :ref:`Geometry <geometry>` documentation.  The ``regenerate`` method does as its name suggests and regenerates all the data of the object.  As methods are added to the object (using ``add_method``) they are noted in a list, which is then referenced at regeneration time so all custom added methods are invoked. 

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Phases
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The **Phases** class contains methods for estimating or predicting the thermo-physical properties of fluids, such as viscosity or density, and this class can also be used to calculate solid properties.  (The name 'phase' was chosen over 'fluid' to be more general).

**Phase** objects are 'built' by the user to contain the specific methods that are to be used to calculate the phase properties.  For instance, a **Phase** object can calculate viscosity assuming a constant value, or Reynolds equation.  The user can use the methods supplied with OpenPNM, or add their own.  

Typically there will be multiple **PHase** objects defined for each simulation, since most models will have at least an invading fluid and a defending fluid.  There can be an unlimited number of phases associated with a **Network**.  

This module works identically to the `Geometry`_ module.    **Phase** objects are essentially empty when initialized, and the user adds the desired methods.  For example, the molar density of a fluid is a fundamental property that is required for many other calculations.  OpenPNM comes with several models that can calculate the molar density of a fluid.  These stored in a submodule called **molar_density**, and options include ``ideal_gas`` and ``real_gas``.  It is envisioned that each user may add their own methods for generating seeds (or ay other calculation) to this submodule (or any other), so the process of adding methods to the object was designed to be as flexible and customizable as possible.  For more information on this process see the :ref:`fluids <Fluids>` documentation.  The ``regenerate`` method does as its name suggests and regenerates all the data of the object.  As methods are added to the object (using ``add_method``) they are noted in a list, which is then referenced at regeneration time so all custom added methods are invoked. 

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Physics
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
asdf

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Algorithms
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
asdf

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Visualization
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
asdf

