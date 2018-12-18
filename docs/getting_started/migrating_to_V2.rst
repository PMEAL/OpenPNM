.. _migrating_to_V2:

================================================================================
Migrating from V1 to V2
================================================================================

--------------------------------------------------------------------------------
Renaming Modules, and Reorganizing
--------------------------------------------------------------------------------

One of the main objectives of V2 was to 'reorganize' the package.  When we started building OpenPNM in 2012 we were new to Python and made some errors in the layout of the code that we fixed in V2.  Specifically, V2 fully adopts the pep8 naming conventions, which suggests that modules be named with lowercase while classes be uppercase.

.. code-block:: python

    # Old way
    import OpenPNM as op
    pn = op.Network.Cubic(shape=[10, 10, 10])

    # New way
    import openpnm as op
    pn = op.network.Cubic(shape=[10, 10, 10])

Note that the package itself was renamed to lowercase, as well as all the module names.  The tables below give an overview of the package comparing V1 and V2:

The package itself was renamed using all lowercase letters:

+--------------------------------------+---------------------------------------+
| Version 1                            | Version 2                             |
+======================================+=======================================+
| OpenPNM                              | openpnm                               |
+--------------------------------------+---------------------------------------+

All the submodules were renamed using lowercase letters:

+--------------------------------------+---------------------------------------+
| Version 1                            | Version 2                             |
+======================================+=======================================+
| Network                              | network                               |
+--------------------------------------+---------------------------------------+
| Phases                               | phases                                |
+--------------------------------------+---------------------------------------+
| Geometry                             | geometry                              |
+--------------------------------------+---------------------------------------+
| Physics                              | physics                               |
+--------------------------------------+---------------------------------------+
| Algorithms                           | algorithms                            |
+--------------------------------------+---------------------------------------+
| Utilities                            | utils                                 |
+--------------------------------------+---------------------------------------+

The terms Base and Core were switched; the Base module is now called core and the Core class is now Base.  This is more consistent with other packages, where all the 'core' classes are stored in a core module, and the 'Base' class is the main class from which all other classes inherit (e.g. Cubic <-- GenericNetwork <-- Base):

+--------------------------------------+---------------------------------------+
| Version 1                            | Version 2                             |
+======================================+=======================================+
| Base                                 | core                                  |
+--------------------------------------+---------------------------------------+
| Base.Core                            | core.Base                             |
+--------------------------------------+---------------------------------------+

Pore-scale models are now stored in a top-level module to improve findability:

+--------------------------------------+---------------------------------------+
| Version 1                            | Version 2                             |
+======================================+=======================================+
| Phases.models                        | models.phases                         |
+--------------------------------------+---------------------------------------+
| Geometry.models                      | models.geometry                       |
+--------------------------------------+---------------------------------------+
| Physics.models                       | models.physics                        |
+--------------------------------------+---------------------------------------+

A new module was added called topotools to house all functions pertaining to manipulating and mangling the topology of the networks.  All functions formerly found under models and/or tools have been moved:

+--------------------------------------+---------------------------------------+
| Version 1                            | Version 2                             |
+======================================+=======================================+
| Network.models                       | *removed*                             |
+--------------------------------------+---------------------------------------+
| Network.tools                        | topotools                             |
+--------------------------------------+---------------------------------------+

Finally, the Postprocessing module was removed and the few functions it had worth keeping have been relocated:

+--------------------------------------+---------------------------------------+
| Version 1                            | Version 2                             |
+======================================+=======================================+
| Postprocessing                       | *removed*                             |
+--------------------------------------+---------------------------------------+

--------------------------------------------------------------------------------
Cleaning up Classes
--------------------------------------------------------------------------------

In V2, we have made an effort to reduce the number of methods on each class, particularly GenericNetwork. Many of the functions that were previously found on the various networks (trim, extend) are now stored in topotools, and instead of being class methods are simple functions that accept a network as an argument an operate on the network 'in-place' meaning they do not return a network, but the received network is changed.

For instance:

.. code-block:: python

    # Old way:
    pn = op.Network.Cubic(shape=[10, 10, 10])
    pn.trim(pores=[1, 4, 7])

    # New way:
    pn = op.network.Cubic(shape=[10, 10, 10])
    op.topotools.trim(network=pn, pores=[1, 4, 7])

--------------------------------------------------------------------------------
Instantiating Physics Objects
--------------------------------------------------------------------------------

In version 2 we have moved to a more strict assignment of Physics objects to different regions, so instead of assigning them to specific pores and throats (as is still done for Geometry objects), each Physics is associated 1-to-1 with a Geometry.  For example:

.. code-block:: python

    # Old way
    import OpenPNM as op
    pn = op.Network.Cubic(shape=[10, 10, 10])
    geom = op.Geometry.GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)
    water = op.Phases.Water(network=pn)
    phys = op.Physics.GenericPhysics(network=pn, phase=water, pores=pn.Ps, throats=pn.Ts)

    # New way:
    import openpnm as op
    pn = op.network.Cubic(shape=[10, 10, 10])
    geom = op.geometry.GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)
    water = op.phases.Water(network=pn)
    phys = op.physics.GenericPhysics(network=pn, phase=water, geometry=geom)

The changes to the Physics instantiation were motivated by a new way organize the simulations called "The Grid".  The grid is a 2D table of objects where each column header is a Phase, each row label is a Geometry, and each row-column intersection contains the corresponding Physics.  By forcing a 1-to-1 associated with Geometry and Physics, we maintain the grid structure which is helpful for object look-ups.  The grid can be seen using the new Project object as follows:

.. code-block:: python

    import openpnm as op
    pn = op.network.Cubic(shape=[10, 10, 10])
    geom = op.geometry.GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)
    water = op.phases.Water(network=pn)
    phys = op.physics.GenericPhysics(network=pn, phase=water, geometry=geom)
    proj = pn.project  # (or geom.project, or water.project or phys.project)
    print(proj.grid)

Which will produce the following:

::

    ――――――――――――――――――――――――――――――――
    |    net_01    |   phase_01    |
    ――――――――――――――――――――――――――――――――
    |    geo_01    |    phys_01    |
    ――――――――――――――――――――――――――――――――

--------------------------------------------------------------------------------
Object Lookups
--------------------------------------------------------------------------------

Object look-ups have changed substantially.  In V1 every object possessed a handle to it's associated objects, such that ``water.network`` would actually contain a handle directly to the network object. This association was causing problems such as memory leaks, large objects with circular references, and complicated object deletion.  In V2 we have moved away from this completely, and now all object look-ups are done with a "top down" approach.  This means that you must ask the Project object to find the associated object.  For example:

.. code-block:: python

    # Old way
    phys.phases()

    # New way
    proj = phys.project
    proj.find_phase(physics=phys)

The actual looking up is done by checking the labels on each object.  When an object such as ``phys`` was is created, it is given a name, such as ``phys_01``, and the Phase with which it is associated is updated to have the labels ``'pore.phys_01'`` and ``'throat.phys_01'`` with ``True`` values where the ``phys`` object was assigned.  This 'top-down' object look up will search each Phase object until it finds the label corresponding to ``phys``, and once it does, the associated Phase is found.  As similar approach is done to find which geometry is associated with which physics and so forth.  

--------------------------------------------------------------------------------
Dealing with Pore-Scale Models
--------------------------------------------------------------------------------

Pore-scale models are one of the more confusing aspects for new users of OpenPNM, so we have attempted to simplify them in several ways.  Each object still has a ``models`` attribute which is a dictionary containing all the model definitions, but this dictionary, called the ModelsDict, is much simpler to work with.

1.  Models are now run in order of dependencies, so if a user adds a model for 'pore.volume' before adding 'pore.diameter', OpenPNM will detect this and run the 'pore.diameter' model first so the values are available to subsequent models.  In the event that OpenPNM cannot run models in the correct order (e.g. a needed model is missing), it will report a warning, and can be remedied by adding the missing model and calling ``regnerate_models``.

2.  ``add_model`` and ``regnerate_models`` are now methods of the Base class rather than being methods of the models dict.  This was done to increase their visibility and make them easier to call.  For example:

.. code-block:: python

    # Old way
    geom.models.add(<args go here>)
    geom.models.regnerate()

    # New way
    geom.add_model(<args go here>)
    geom.regnerate_models()

3.  Models themselves are now much simpler and more flexible.  Any given model can be assigned to any object. The aim where was to allow, for instance, geometrical models like 'pore.diameter' to be assigned directly to a Network object.  The point of Geometry objects is to allow complex, heterogeneous materials (e.g. layers with different pore sizes), but this is not always needed, so users can bypass this feature.  Having no Geometry objects precludes the use of Physics objects (see previous section on instantiating Physics objects), which means that pore-scale physics models (e.g. capillary pressure) must be assigned to Phase objects, which is also possible under this new flexible scheme.  In all pore scale models, the object to which the model is attached is referred to as the ``target``.

--------------------------------------------------------------------------------
Other Breaking Changes
--------------------------------------------------------------------------------

1. Neumann and Dirichlet boundary conditions have been renamed to Rate and Value to be more descriptive.  Neumann in particular was a problem since we were not actually specifying the gradient in the BC, but the gradient multiplied by the conductance, hence the move to Rate.  Dirichlet was renamed to Value was for consistency with Rate.

2. The behavior of the various lookup methods (pores, find_neighbor_pores) now all use the same ``mode`` keywords. The new keyswords are taken for logic theory and include 'and', 'or', 'xor', 'xnor'.  The functions also accept synonyms from set theory ('intersection', 'union', 'exclusive_or', while xnor has no equivalent).

.. note:: The meaning of intersection has changed

  **One very important change** is the meaning of ``mode='intersection'`` in the ``find_neighbor_pores`` method.  In version 1 this was incorrectly being used to find neighbors that shared *more than one* input pore, while in version 2 this means pores that share *all the* input pores.  The mode 'xnor' replaces the old 'intersection'.  This change needed to be made, but is problematic since 'intersection' is still an accepted mode but returns a different result

3. Every object now has a ``settings`` attribute, which is a dictionary of key:value pairs.  This was added so that we could stop storing various flags and values as attributes (sometimes hidden) on objects, and keep all such information in one place.  All Algorithm objects in particular use this settings dictionary, and the ``setup`` method essentially just passes any received arguments onto the ``settings`` dict.
