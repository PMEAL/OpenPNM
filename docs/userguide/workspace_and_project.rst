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

Note: Each Project can contain one and only one Network.  Since every other object must be associated with a single network, so it follows that there is only one network per project.

--------------------------------------------------------------------------------
Using Workspace and Project Objects
--------------------------------------------------------------------------------

The Workspace and the Project work together.  The Workspace is the highest level object which tracks all of the open Projects.  Projects provide a lower level of oversight, with each Project keeping of track of the specific OpenPNM objects for each simulation.

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
Removing and Moving Objects in a Project
................................................................................

Removing an object from a Project can be done with the ``purge_object`` method
or with the ``remove`` method of the ``list`` class (which has been subclassed
to be a wrapper for ``purge_object``).

.. code-block:: python

    >>> geo1 in proj1
    True
    >>> proj1.purge_object(geo1)
    >>> geo1 in proj1
    False

It's also quite easy to copy objects between projects.  For e










.
