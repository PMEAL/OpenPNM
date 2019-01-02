.. _quick_start:

================================================================================
Quick Start
================================================================================

.. contents:: Page Contents
    :depth: 3

The following is meant to give a very quick overview of how OpenPNM works.

--------------------------------------------------------------------------------
Creating a Network
--------------------------------------------------------------------------------

The first step in an OpenPNM project is to create a network.

.. code-block:: python

    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[10, 10, 10], spacing=0.0001)

The resulting network can be quickly visualized using ``plot_coordinates`` and
``plot_connections`` in the :ref:`topotools_index` module to get the following:

.. image:: /../docs/static/images/quick_start_cubic_network.png
    :width: 800px
    :align: center

The ``network`` module has a number of network types to chose from, and they can be found on the :ref:`network_index` page. You can also import networks from various outside sources, include networks that have been extracted from tomographic images using our `PoreSpy Package <http:\\porespy.com>`_.

--------------------------------------------------------------------------------
Adding Geometrical Properties
--------------------------------------------------------------------------------

The Network only contains spatial information (pore coordinates) and topological information (throat connections), so it is necessary to add geometrical information by creating a Geometry object:

.. code-block:: python

    >>> Ps = pn.pores('all')
    >>> Ts = pn.throats('all')
    >>> geo = op.geometry.StickAndBall(network=pn, pores=Ps, throats=Ts)

In this case the `StickAndBall` class was used, which has preset pore-scale models that calculate properties such as pore diameters and throat lengths, based on the assumption that the pores are spherical and the throats are cylinders.

.. note::

    (1) The Network (``pn``) was passed to the Geometry class as an argument, so that ``geo`` knows which network it's associated with.

    (2) the Geometry was assigned to specified pores (``Ps``) and throats (``Ts``), in this case it was 'all' of them but it's possible to use several different Geometry objects for different subsets of the domain.

--------------------------------------------------------------------------------
Creating Phases
--------------------------------------------------------------------------------

Phases must be created to calculate the thermophysical properties of the fluids (and solids) used in the simulations:

.. code-block:: python

    >>> hg = op.phases.Mercury(network=pn)
    >>> h2o = op.phases.Water(network=pn)

OpenPNM includes a few common phases, including :ref:`mercury_api`, :ref:`air_api` and :ref:`water_api`, but also a set of pore-scale models for calculating properties of different phases.

.. note::

    Phase objects are associated with a Network, but they are not assigned to specific pores and throats.  This is because phases can exist anywhere and everywhere in the domain, and can move around.

--------------------------------------------------------------------------------
Assigning Pore-Scale Physics Models
--------------------------------------------------------------------------------

.. code-block:: python

    >>> phys_hg = op.physics.GenericPhysics(network=pn, phase=hg, geometry=geo)
    >>> phys_h2o = op.physics.GenericPhysics(network=pn, phase=h2o, geometry=geo)

The ``GenericPhysics`` class was used, which has NO pore-scale models attached.  We will add this manually in the next step.

.. note::

    (1) The Network must be given as an argument so that the Physics knows which network it's associated with

    (2) One Physics object is required for each Phase, since physics models require thermophysical properties.  For example, the Hagan-Poisseuille equation requires the viscosity of the phase.

    (3) Each Physics object is also associated with a Geometry.  The reason for this is to assign different pore-scale physics models to different regions.  In this case both are associated with ``geo`` since there is only one Geometry in the whole domain.

................................................................................
Assigning Pore-Scale Models
................................................................................

We must assign models to each of our Physics.  The ``hg`` phase will be used to simulate a Porosimetry experiment, so it needs a capillary pressure model, and ``h2o`` will be used in a permeability simulation so we must define a hydraulic conductance.   The following shows how to fetch models from the ``models`` library, attach them to the target object, and specify the model parameters:

.. code-block:: python

    >>> model = op.models.physics.capillary_pressure.washburn
    >>> phys_hg.add_model(propname='throat.entry_pressure',
    ...                   model=model,
    ...                   contact_angle='pore.contact_angle',
    ...                   surface_tension='pore.surface_tension')
    >>> model = op.models.physics.hydraulic_conductance.hagen_poiseuille
    >>> phys_h2o.add_model(propname='throat.hydraulic_conductance',
    ...                    model=model,
    ...                    pore_viscosity='pore.viscosity',
    ...                    pore_area='pore.area',
    ...                    throat_area='throat.area',
    ...                    conduit_lengths='throat.conduit_lengths')

--------------------------------------------------------------------------------
Performing Some Simulations
--------------------------------------------------------------------------------

We are now ready to conduct some simulations.  The most important step in validating a pore network model is to ensure that it reproduces experimentally measured porosimetry curves and absolute permeability.

................................................................................
Simulating Mercury Porosimetry
................................................................................

.. code-block:: python

    >>> mip = op.algorithms.Porosimetry(network=pn)
    >>> mip.setup(phase=hg)
    >>> mip.set_inlets(pn.pores(['left', 'right', 'top', 'bottom', 'front',
    ...                          'back']))
    >>> mip.run(points=25)


Which can be visualized using the ``plot_intrusion_curve`` method of the Porosimetry class:

.. image:: /../docs/static/images/quick_start_drainage_curve.png
    :width: 800px
    :align: center

................................................................................
Calculating Network Permeability
................................................................................

Similarly for the permeability calculation:

.. code-block:: python

    >>> perm = op.algorithms.StokesFlow(network=pn)
    >>> perm.setup(phase=h2o)
    >>> perm.set_value_BC(pores=pn.pores('left'), values=1)
    >>> perm.set_value_BC(pores=pn.pores('right'), values=0)
    >>> # perm.run()

The above code solves for the pressure in each pore and stores the result as ``perm['pore.pressure']``.  To find the permeability of the network, there is a ``calc_eff_permeability`` method on the StokeFlow class.  Start by telling the algorithm the area and length of the domain (unfortunately there is no sure way to get these accurately, though values will be guesstimated if not provided)

    >>> perm.domain_area = (10*0.0001)**2
    >>> perm.domain_length = (10*0.0001)

Finally to run the algorithm and calculate network permeability use:

.. code-block:: python

    perm.run()
    K = perm.calc_effective_permeability()

.. note::

    (1) The ``calc_eff_permeability`` method finds K by inverting Darcy's law, and looking up all the necessary information (pressure drop, viscosity) from the objects which the algorithm is associated.

    (2) If the domain area and length are not given, an attempt is made to estimate them but it's more accurate to provide it.

--------------------------------------------------------------------------------
Saving Project and Exporting to Paraview
--------------------------------------------------------------------------------

Now that the simulation is finished, it can be saved to a ``.pnm`` file for future use.  OpenPNM has two levels of *management*: the Workspace and the Project.  Each Project contains a single network and its associated object (all the code in this guide are a single Project).  The Workspace contains all the active Projects.  You can save the entire Workspace including all active Projects, or you can save a single Project.

Each object has a ``project`` attribute which returns a handle to the Project to which it belongs.

.. code-block:: python

    >>> proj = pn.project  # Retrieve the project handle

The Project object offers several useful tools, including the ability to ``export_data`` to various formats, including VTK for viewing in `Paraview <http://www.paraview.org>`_.  Using Paraview provides much better visualization than the ``plot_connections`` and ``plot_coordinates`` used above:

.. image:: http://i.imgur.com/GbUNy0b.png
   :width: 500 px
   :align: center
