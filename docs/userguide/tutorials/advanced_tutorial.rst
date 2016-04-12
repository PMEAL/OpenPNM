.. _advanced_tutorial:

###############################################################################
Tutorial 3 of 3: Advanced Topics and Usage
###############################################################################

In this final tutorial on the use of OpenPNM, simulate pressure drive gas flow through a porous medium that is partially invaded by a non-wetting fluid.  This is known as relative permeability and is one of the standard simulations performed with pore networks.  This

**Learning Outcomes**

#. Experiment with adding boundary pores
#. Manipulate Network Topology
#. Explore the ModelsDict design, including copying models between objects, and changing model parameters
#. Combine multiple algorithms to predict relative permeability
#. Write a custom pore-scale model
#. Use the workspace manager to save and load, clone and purge simulations

===============================================================================
Build Network Topology
===============================================================================

For the present tutorial, we'll keep the topology simple to help keep the focus on other aspects of OpenPNM.

.. code-block:: python

    >>> import OpenPNM as op
    >>> pn = op.Network.Cubic(shape=[10, 10, 10], spacing=0.00006)

-------------------------------------------------------------------------------
Adding Boundary Pores
-------------------------------------------------------------------------------

When performing transport simulations it is often useful to have 'boundary' pores attached to the surface(s) of the network where boundary conditions can be applied.  The **Cubic** class has two methods for doing this: ``add_boundaries`` and ``add_boundary_pores``.  The first method adds boundaries to ALL six faces of the network and offsets them from the network by 1/2 of the value provided as the network ``spacing``.  The second method provides total control over which boundaries are created and where they are positioned, but it more cumbersome to use.  Let's explore these options:

.. code-block:: python

    >>> Np = pn.Np  # Should be 1000
    >>> Nt = pn.Nt  # Should be 2700
    >>> pn.add_boundaries()
    >>> Np2 = pn.Np  # Should be 1600 (10x10 pores per face)
    >>> Nt2 = pn.Nt  # Should be 3300 (1 throat connecting each new pore)

Let's remove all these newly created boundary pores.  When they are created these pores are all automatically labeled with a label such as ``'top_boundary'``, so we can select all boundary pores using:

.. code-block:: python

    >>> Ps = pn.pores('*boundary')  # Using the * wildcard

We can then ``trim`` these pores from the network using:

.. code-block:: python

    >>> pn.trim(pores=Ps)
    >>> Np = pn.Np  # Should be 1000 now
    >>> Nt = pn.Nt  # Should be 2700

Note that all throats connecting to the trimmed pores were automatically removed since OpenPNM does not allow 'dangling' or 'headless' throats.

Now that ``pn`` is back to its original size, let's explore the second approach to apply boundary pores.

.. code-block:: python

    >>> Ps = pn.pores('top')  # Select pores on top of network
    >>> pn.add_boundary_pores(pores=Ps, offset=[0, 0, 0.00003],
    ...                       apply_label='top_boundary')
    >>> Ps = pn.pores('bottom')  # Select pores on bottom of network
    >>> pn.add_boundary_pores(pores=Ps, offset=[0, 0, -0.00003],
    ...                       apply_label='bottom_boundary')
    >>> Np = pn.Np  # Should be 1200 (10x10 for two faces)
    >>> Nt = pn.Nt  # Should be 2900

This approach requires more typing than the ``add_boundaries`` method, but allows for much finer control over how boundaries are created.

===============================================================================
Define Geometry Objects
===============================================================================

Since we've added boundary pores to the network we need to the treat them a little bit differently.  Specifically, they should have no volume or length (as they are not physically representative of real pores).  To do this, we create two separate **Geometry** objects, one for internal pores and one for the boundaries:

.. code-block:: python

    >>> Ps = pn.pores('*boundary', mode='not')
    >>> Ts = pn.find_neighbor_throats(pores=Ps, mode='intersection')
    >>> geom = op.Geometry.Stick_and_Ball(network=pn, pores=Ps, throats=Ts)
    >>> Ps = pn.pores('*boundary')
    >>> Ts = pn.find_neighbor_throats(pores=Ps)
    >>> boun = op.Geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)

The **Stick_and_Ball** class is preloaded with the pore-scale models to calculate all the necessary size information (diameter, lengths, etc).  The **GenericGeometry** class used for the boundary pores and throats is empty and requires our input.  Since boundary pores are fictitious we want them to have suitable properties:

.. code-block:: python

    >>> boun['pore.diameter'] = 0
    >>> boun['pore.volume'] = 0

Boundary throats act as the link between the internal pores and the 'outside', and hence should be considered as real throats.  For this, we will add some pore-scale models:

.. code-block:: python

    >>> boun.models.add(propname='throat.length',
    ...                 model=op.Geometry.models.throat_length.straight)
    >>> boun.models.add(propname='throat.diameter',
    ...                 model=op.Geometry.models.throat_misc.neighbor,
    ...                 pore_prop='pore.diameter')  # More on this model below
    >>> boun.models.add(propname='throat.area',
    ...                 model=op.Geometry.models.throat_area.cylindrical)
    >>> boun.models.add(propname='throat.volume',
    ...                 model=op.Geometry.models.throat_volume.cylinder)

These models are required for the Hagan-Poiseiulle model. Most of them are straight-forward geometry calculations, except for the model used for ``'throat.diameter'``.  In this case the model looks into the neighbor pores, retrieves the two ``'pore.diameter'`` and uses the ``'max'`` value.  Because we set the boundary pores to have 0 diameter, this will naturally find result in the throat being assigned the diameter of the internal pore.
