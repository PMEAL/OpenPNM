.. _models_guide:

================================================================================
Pore Scale Models
================================================================================

.. contents:: Page Contents
    :depth: 3

Models are one of the most important aspects of OpenPNM, as they allow the user to specify a 'model' for calculating properties (e.g. 'pore.volume'), rather than just entering numerical values (e.g. ``geometry_object['pore.volume'] = array`` ).

Models offer several vital functions:

1.  The ability to save recipes
2.  The ability to regenerate all values when one changes
3.  Creating custom behavior

--------------------------------------------------------------------------------
The Models Module: A Library of Prewritten Models
--------------------------------------------------------------------------------

OpenPNM comes with a wide assortment of models for calculating all sorts of things such as geometrical properties, thermophysical fluid properties, and physics parameters.  These are stored under the ``models`` attribute and categorized by the type of properties they produce.  For instance, *'openpnm.models.geometry.pore_diameter'* contains several methods for calculating pore diameters.

Most `Base objects <base_api>`_ (`Geometry <generic_geometry_api>`_, `Phases <generic_phase_api>`, `Physics <generic_physics_api>`_) have a ``models`` attribute which upon instantiation of the object is filled with an empty `ModelsDict <modelsdict_api>`_ object.  The **ModelsDict** object is designed to store and interact with all *models* on the **Base** object.  The **ModelsDict** is a subclass of the Python *dictionary* type, with several features added for dealing specifically with *models*.  Each *model* and its associated arguments are stored under a single dictionary key, ** under the specified *'propname'*.  When a model is run the values it produces are automatically stored in the **Core** object's dictionary under the same specified *'propname'*.

--------------------------------------------------------------------------------
Adding Models
--------------------------------------------------------------------------------

Adding a model to an object is done as follows:

(1) A handle to the desired model is retrieved, either from the included OpenPNM model libraries, or from a file containing the users custom models.
(2) The model is attached to the target object using ``add_model``.

This process is demonstrated by adding a random pore seed model to a **Geometry** object:

.. code-block:: Python

    >>> import openpnm as op
    >>> pn = op.network.Cubic([5, 5, 5])
    >>> geom = op.geometry.GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)
    >>> mod = op.models.geometry.pore_size.random  # Get a handle to the desired model
    >>> geom.add_model(propname='pore.seed', model=mod, seed=0)

The *'propname'* and *'model'* arguments are required by the ``add_model`` method, and all other arguments such *'seed'* are passed on the model (In this case it specifies the initialization value for the random number generator).

One can inspect all of the models stored on a given **Base** object by typing ```print(geom.models)``` at the command line:

.. code-block::  Python

    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    #   Property Name             Parameter                 Value
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    1   pore.seed                 model:                    random
                                  seed:                     0
                                  num_range:                [0, 1]
                                  regeneration mode:        normal
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――

--------------------------------------------------------------------------------
Regenerating Models
--------------------------------------------------------------------------------

In order to recalculate the data, the models stored in the **ModelsDict** dictionary must be rerun.  This is accomplished with the ``regenerate_models`` method.  This method takes an optional list of *'propnames'* that should be regenerated.  Models are regenerated in an order determined automatically by OpenPNM to ensure that all dependent properties are calculated first (e.g. 'pore.diameter' is recalculated before 'pore.volume' which is a function of diameter).
