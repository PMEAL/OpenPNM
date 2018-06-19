.. _workspace:

pyreverse -AS -o png -p openpnm "C:\Users\Jeff\Dropbox\Flash Sync\Code\Git\OpenPNM\openpnm\core"

##############################################################################
Workspace Manager and Project Container
##############################################################################

==============================================================================
The Workspace
==============================================================================
OpenPNM includes a **Workspace** Manager object that performs many of the functions found in the *menu bar* of a typical application's GUI, such as saving and loading sessions.

The **Workspace** class is a `Singleton <https://en.wikipedia.org/wiki/Singleton_pattern>`_ in Object Oriented Programming jargon, meaning that only ONE instance can exist at any given time.  In other words, each time a Singleton is instantiated it returns the already existing object if one already exists.  This behavior is handy since it means you can instantiate the **Workspace** at any time, from anywhere in your workflow, and you'll have access to the one and only Workspace object.

==============================================================================
The Project
==============================================================================
Another workflow management tool is the Project object.  A Project defines a single simulation, which would contain a network, some geometry, phase and physics objects, and any simulations.  The main roles of the Project is to group related objects together and to provide additional administrative tools.

The Project object is a Python ``list`` that has been subclassed to have many additional methods and functions, purging an object or finding a list of all phase object in the current project.

Note: Each Project can contain one and only one Network.  Since every other object must be associated with a single network, so it follows that there is only one network per project.

==============================================================================
Using Workspace and Project Objects
==============================================================================

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
    >>> pn = op.network.Cubic(shape=[3, 3, 3], project=proj1)
    >>> pn in proj1
    True
