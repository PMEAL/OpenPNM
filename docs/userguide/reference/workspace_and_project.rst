.. _workspace:

pyreverse -AS -o png -p openpnm "C:\Users\Jeff\Dropbox\Flash Sync\Code\Git\OpenPNM\openpnm\core"

##############################################################################
Workspace Manager
##############################################################################
OpenPNM includes a **Workspace** class that performs many of the functions found in the *menu bar* of a typical application's GUI, such as saving and loading sessions.

The **Workspace** class is a `Singleton <https://en.wikipedia.org/wiki/Singleton_pattern>`_ in Object Oriented Programming Jargon, meaning that only ONE instance can exist at any given time, and moreover, each time a Singleton is instantiated it returns the already existing object if one already exists.  This behavior is handy since it means you can instantiate the **Workspace** at any time, from anywhere in your workflow, and you'll have access to the one and only object.

All OpenPNM objects register themselves with this single **Workspace** when they are created, so you can access any existing object via the workspace.  Like so many custom classes in Python, the **Workspace** is a *dictionary* and each OpenPNM object is stored by name.  For example:

.. code-block:: python

    >>> import openpnm as op
    >>> mgr = op.core.Workspace()
    >>> mgr.clear()  # Clear workspace of any pre-existing objects
    >>> mgr.keys()
    dict_keys([])
    >>> pn = op.network.Cubic(shape=[5, 5, 5], name='foo')
    >>> list(mgr.keys())
    ['foo']
    >>> pn2 = op.network.Cubic(shape=[5, 5, 5], name='bar')
    >>> sorted(list(mgr.keys()))
    ['bar', 'foo']
    >>> pn is mgr['foo']  # The object stored as 'foo' actually pn
    True

The **Workspace** object also tracks the relationships between the OpenPNM objects, and the ``print`` function outputs a formatted list of each simulation's structure:

.. code-block:: python

    >>> geom = op.Geometry.GenericGeometry(network=pn, name='geom_on_foo')
    >>> geom2 = op.Geometry.GenericGeometry(network=pn2, name='geom_on_bar')

.. code-block:: python

    print(mgr)
    ------------------------------------------------------------
    Object:         Name                 (Class)
    ------------------------------------------------------------
    Network:        bar                  (Cubic)
    ++ Geometry:    geom_on_bar          (GenericGeometry)
    ------------------------------------------------------------
    Object:         Name                 (Class)
    ------------------------------------------------------------
    Network:        foo                  (Cubic)
    ++ Geometry:    geom_on_foo          (GenericGeometry)

A list of all the methods available on the **Workspace** object can be obtained with:

.. code-block:: python

    >>> methods = [item for item in dir(mgr) if not item.startswith('_')]

This list contains the usual *dictionary* methods, plus an assortment of others for specifically managing OpenPNM objects.  The sections below will discuss each in more detail.

===============================================================================
Save and Load the Entire Workspace
===============================================================================

If for some reason you want to save ALL of your current objects, then you can use ``mgr.save``.  This will create a ``.pnm`` file in your current working directory with all objects faithfully saved (including data, labels, models, model arguments, relationships, everything.).  This is accomplished by the `Python Pickle library <https://docs.python.org/3/library/pickle.html>`_.  If no name is supplied then a name will automatically be generated using the current date and time.

.. autosummary::

    ~OpenPNM.Base.Controller.save

Loading a simulation that was previously saved to a file is equally straightforward, using ``mgr.load``.  Of course, a file name must be given.  When loading a simulation from a file, any objects stored in the **Workspace** will be deleted.

.. autosummary::

    ~OpenPNM.Base.Controller.load

.. note::

    Saving and loading entire workspaces is really only useful for simulations that take a long time to generate, such as **Delaunay** networks with **Voronoi** geometry.  For fast simulations it is just as easy to save the *'.py'* script, then to recreate a whole new simulation on demand.

.. warning::

    **Algorithm** objects are not automatically registered with the **Workspace** when they are created.  This is because in some cases algorithms are instantiated inside a *for-loop* which would quickly bloat the size of the *'.pnm'* file.  This may change in a future version.


===============================================================================
Save and Load Individual Simulations
===============================================================================

Instead of saving the entire workspace it is also possible to save individual simulations.  For instance, if multiple networks have been defined but only one of them is of interest, then that **Network** along with all the **Geometry**, **Phase**, and **Physics** objects which were associated with it can be saved using ``mgr.save_simulation``.

.. code-block:: python

    >>> pn1 = op.Network.Cubic(shape=[10, 10, 10])
    >>> geo = op.Geometry.GenericGeometry(network=pn1, pores=pn1.Ps, throats=pn1.Ts)
    >>> pn2 = op.Network.Cubic(shape=[10, 10, 10])
    >>> geo2 = op.Geometry.GenericGeometry(network=pn2, pores=pn2.Ps, throats=pn2.Ts)
    >>> air = op.Phases.Air(network=pn1)
    >>> water = op.Phases.Water(network=pn2)
    >>> mgr.save_simulation(network=pn1, filename='first_network.net')
    >>> mgr.save_simulation(network=pn2, filename='second_network.net')

The above lines create two files in the current working directory called *'first_network.net'* and *'second_network.net'* which contain pn1 and pn2 respectively, along with all objects (ie. **Geometry** and **Phase**) associated with each.

If we now ``clear`` the **Workspace** object, we can reload each of these simulations:

.. code-block:: python

    >>> mgr.clear()
    >>> mgr.load_simulation('first_network.net')
    >>> mgr.load_simulation('second_network.net')

When loading multiple 'simulations' into the **Workspace** it does not remove any existing simulations (unlike loading a saved workspace *'.pnm'* file).

The ``save_simulation`` and ``load_simulation`` methods are ideal when running large batches of calculations and you want to save the numerical results for later analysis.

.. warning::

    **Algorithm** objects are not automatically registered with the **Workspace** when they are created.  This is because in some cases algorithms are instantiated inside a *for-loop* which would quickly bloat the size of the *'.net'* file.  This may change in a future version.

===============================================================================
Import and Export Data
===============================================================================

The **Workspace** manager has methods to ``import_data`` and ``export_data``.  These are wrapper methods for the actual methods found in ``Utilities.IO``.  These wrapper or helper methods accept several arguments that control which type of file is imported or exported.  The actual import and export is explained fully in :ref:`data_IO`.

===============================================================================
Object Lookup
===============================================================================

Each OpenPNM Core object that is created is either given or assigned a ``name``.  This name is used as the dictionary key when the object is saved on the **Workspace** manager, as outlined above.  In addition to looking up objects by name, it is also possible to look them up by type using ``networks``, ``geometry``, ``physics``, and ``phases``.  At present ``algorithms`` is offered but does not return any objects since **Algorithms** are not registered.  Each of these methods returns a *list* of objects of the specified type.  The objects in the list can be assigned to variables on the command line for easy access:

.. code-block:: python

    >>> pn = mgr.networks()[0]

===============================================================================
Object Manipulation (Purging, Cloning, etc)
===============================================================================

::

    Documentation not finished yet
