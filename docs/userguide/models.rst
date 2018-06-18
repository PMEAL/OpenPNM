.. Modelling multiphase transport in fuel cells: The power of pore-scale approaches:

===============================================================================
Pore Scale Models
===============================================================================
Models are one of the most important aspects of OpenPNM, as they allow the user to specify a 'model' for calculating 'pore.volume', rather than just entering numerical values into a ``geometry_object['pore.volume']`` array for instance.  It is mainly through customized models that users can tailor OpenPNM to a specific situation, though OpenPNM includes a variety of pre-written models.  These are stored under each Module in a folder called 'models'.  For instance, *'Geometry.models.pore_diameter'* contains several methods for calculating pore diameters.  For an example on creating custom models see :ref:`Customizing OpenPNM<customizing>`.

There are two special classes defined in OpenPNM for working with *models*: the **ModelsDict** and the **ModelWrapper**.  Each of these are explained in the following sections.  A schematic of their interrelation with each other and the **Core** object to which they are attached is given below:

.. image:: http://i.imgur.com/WhOY6o7.png

Each **Core** object has a ``models`` attribute which upon instantiation of the **Core** object is filled with an empty **ModelsDict** object.  The **ModelsDict** object is designed to store and interact with all *models* on the **Core** object.  The **ModelsDict** is a subclass of the Python *dictionary* type, with several features added for dealing specifically with *models*.  Each *model* and its associated arguments are wrapped in a **ModelWrapper** dictionary, then added to the **ModelsDict** under the specified *'propname'*.  When a model is run the values it produces are automatically stored in the **Core** object's dictionary under the same specified *'propname'*.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ModelsDict: The Collection of All Models
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Adding a model to an object is done as follows:

(1) A handle to the desired model is retrieved, either from the included OpenPNM model libraries, or from a file containing the users custom models.
(2) The model is attached to the target object using ``add_model``.

This process is demonstrated by adding a random pore seed model to a **Geometry** object:

.. code-block:: Python

    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,pores=pn.Ps,throats=pn.Ts)
    >>> mod = OpenPNM.Geometry.models.pore_misc.random  # Get a handle to the desired model
    >>> geom.add_model(propname='pore.seed', model=mod, seed=0)

The *'propname'* and *'model'* arguments are required by the ``add_model`` method, and all other arguments such *'seed'* are passed on the model (In this case it specifies the initialization value for the random number generator).

One can inspect all of the models stored on a given **Core** object by typing ```print(geom.models)``` at the command line:

.. code-block::  Python

    ------------------------------------------------------------
    #     Property Name                  Regeneration Mode
    ------------------------------------------------------------
    1     pore.seed                      normal
    ------------------------------------------------------------

By default, the ``add_model`` method runs the model and places the data in the Core object's dictionary under the given *'propname'*. It also wraps the handle to the model and all the given parameters into a **ModelWrapper** dictionary (described below), then saves it in the **ModelsDict** under the same *'propname'*.  There are several opions for the *regeneration mode* such as *'deferred'* and *'on_demand'*.  Each of these is described in the docstring for the ``add_model`` method.

In order to recalculate the data, the models stored in the **ModelsDict** dictionary must be rerun.  This is accomplished with the ``regenerate`` method.  This method takes an optional list of *'propnames'* that should be regenerated.  Models are regenerated in the order that they were added to the object so some care must be taken to ensure that changes in property values cascade through the object correctly.  The **ModelsDict** class has functions for updating the model order (``reorder``).

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ModelWrapper: The Home for Each Individual Model
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Each of the models that are added to the **ModelsDict** are first wrapped inside a **ModelWrapper**, which is a subclass of Python's *dictionary* type.  Most users will not need to interact with the **ModelWrapper** directly, since most applications only use the ``add_model`` and ``regenerate`` methods of the **Core** object.  There are a few cases where it is required, such as changing any of the initially set arguments.

When the ``add_model`` method of the **Core** object is called, it collects all of the arguments that it receives and stores them in a new instance of a **ModelWrapper** under the appropriate argument name (i.e. *'seed'* is stored under ``ModelWrapper['seed']``).  It is possible to alter any of these values directly, and these changes will be permanent.  One can inspect all the arguments and their current values stored in a **ModelWrapper* by entering ```print(geom.models['pore.seed'])``` at the command line.

.. code-block:: Python

    ------------------------------------------------------------
    OpenPNM.Geometry.models.pore_misc.random
    ------------------------------------------------------------
    Argument Name        Value / (Default)
    ------------------------------------------------------------
    num_range            [0, 1] / ([0, 1])
    regen_mode           normal / (---)
    seed                 0 / (None)
    ------------------------------------------------------------

Once creation of the **ModelWrapper** is complete, its stored in the **ModelsDict** under the specified *'propname'*.  When ``regenerate`` is called, each of the models stored in the **ModelsDict** is run, in order.  When the ``run`` method on the **ModelWrapper** is called, the handle to the model is retrieved from *'models'*, and it is then called with ALL other entries in the **ModelWrapper** sent as keyword arguments.
