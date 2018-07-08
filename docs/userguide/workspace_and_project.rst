.. _workspace:

================================================================================
Workspace Manager and Project Container
================================================================================

.. contents:: Page Contents
    :depth: 3

--------------------------------------------------------------------------------
The Workspace
--------------------------------------------------------------------------------

OpenPNM includes a **Workspace** Manager object that performs many of the functions found in the *menu bar* of a typical application's GUI, such as saving and loading sessions.

The **Workspace** class is a `Singleton <https://en.wikipedia.org/wiki/Singleton_pattern>`_ in Object Oriented Programming jargon, meaning that only ONE instance can exist at any given time.  In other words, each time a Singleton is instantiated it returns the already existing object if one already exists.  This behavior is handy since it means you can instantiate the **Workspace** at any time, from anywhere in your workflow, and you'll have access to the one and only Workspace object.

--------------------------------------------------------------------------------
The Project
--------------------------------------------------------------------------------

Another workflow management tool is the Project object.  A Project defines a single simulation, which would contain a network, some geometry, phase and physics objects, and any simulations.  The main roles of the Project is to group related objects together and to provide additional administrative tools.

The Project object is a Python ``list`` that has been subclassed to have many additional methods and functions, purging an object or finding a list of all phase object in the current project.

.. note::

    Each Project can contain one and only one Network.  Since every other object must be associated with a single network, so it follows that there is only one network per project.  All object initializations can accept either a Network or a Project.


--------------------------------------------------------------------------------
Using Workspace and Project Objects
--------------------------------------------------------------------------------

The Workspace and the Project work together.  The Workspace is the highest level and it tracks all of the open Projects.  Projects provide a lower level of oversight, with each Project keeping of track of the specific OpenPNM objects for each simulation.

The following illustrates the usage of the Workspace and Project objects.

.. code-block:: python

    >>> import openpnm as op
    >>> ws = op.Workspace()
    >>> ws.clear()  # Clear workspace of any pre-existing objects
    >>> ws.keys()
    dict_keys([])
    >>> # Create an empty Project
    >>> proj1 = ws.new_project(name='one')
    >>> # Ensure Project is empty
    >>> proj1
    []
    >>> # Now create a network and associate it with proj1
    >>> pn1 = op.network.Cubic(shape=[3, 3, 3], project=proj1)
    >>> pn1 in proj1
    True

Each Project can only have *one* Network.  If you try to create a second
Network with the same Project, OpenPNM will complain.

Conversely, since each Network *must* be associated with a Project, one is
automatically created if *not* specified:

.. code-block:: python

    >>> pn2 = op.network.Cubic(shape=[3, 3, 3])
    >>> proj2 = pn2.project
    >>> proj1 == proj2
    False

Now that we've successfully created 2 Networks, there will be 2 Projects open
in the Workspace:

.. code-block:: python

    >>> print(ws.keys())
    dict_keys(['one', 'sim_01'])

The second Project was created automatically, and given a default name of
'sim_01'.

When adding other objects, either the Network or the Project can be specified,
which is possible since there is only one Network per Project:

.. code-block:: python

    >>> geo1 = op.geometry.GenericGeometry(network=pn1, pores=pn1.Ps, throats=pn1.Ts, name='geo_01')
    >>> geo2 = op.geometry.GenericGeometry(project=proj2, pores=pn2.Ps, throats=pn2.Ts, name='geo_02')

Projects can fetched from the Workspace by name, and renamed if
desired:

.. code-block:: python

    >>> proj2 = ws['sim_01']
    >>> proj2.name = 'two'
    >>> print(ws.keys())
    dict_keys(['one', 'two'])

................................................................................
Managing Objects Within a Project
................................................................................

The Project object possesses several methods for dealing with the OpenPNM objects it contains.  One of the main uses of the Project is to lookup associated objects.  For instance, given a Physics object (`phys`), you can find which Phase it was associated with using:

.. code-block:: python

    proj = phys.project
    phase = proj.find_phase(physics=phys)

Note that the Project with which each object is associated can be reached from its `project` attribute.

In addition to these lookup methods (others are `find_physics` and `find_geometry`) the project also has the ability to save and load single objects, as well removing objects from the Project.  This latter ability is worth explaining in more detail.  Consider the Grid introduced when explaining :ref:`overall_design`.  When removing an object, it can either result in an empty space on the Grid, or it may be desirable to remove the entire associated row or column, respectively.  The `purge_object` method, therefore, has the ability to remove an isolated object or all of its associated objects.

When an object is purged, not only is it removed from the Project list, but all references to it in other objects in the form of labels (e.g. net['pore.geo_01']) will be removed.

................................................................................
Managing Projects Within a Workspace
................................................................................

The Workspace object possesses methods for working dealing with Project objects, such as saving, loading, closing, and copying them.

Projects are saved with a `.pnm` file format, and this is done using the Python `pickle <https://pymotw.com/3/pickle/index.html>`_ library for serializing objects.  The saved file is actually a `dictionary` that represents a subset of the Workspace, so loading a `.pnm` file is equivalent to unpickling the file then using the Workspace's `update` method to add the contents.  If the contents of the `.pnm` file are a list rather than a dictionary, this also works, so if you manually save a Project as a list (rather then using `save_project`) its still possible to load it using `load_project`.

Another important function of the Workspace is `clone_project`.  As the name suggests, this creates an exact duplicate of a Project and all its objects, but they are unique in memory.  This is useful for creating sub-networks of a master Network to perform small, quick calculations on.  











.
