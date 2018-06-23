.. _quick_start:

================================================================================
Quick Start
================================================================================

.. contents:: Contents of this Page
    :depth: 2


The following is meant to give a very quick overview of how OpenPNM works.

--------------------------------------------------------------------------------
Creating a Network
--------------------------------------------------------------------------------

The first step in an OpenPNM project is to create a network.

.. code-block:: python

    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[10, 10, 10], spacing=0.0001)

The ``network`` module has a number of network types to chose from:

.. currentmodule:: openpnm.network

.. autosummary::
   :nosignatures:

   Cubic
   CubicDual
   CubicTemplate
   Bravais
   Delaunay
   Voronoi
   Gabriel
   DelaunayVoronoiDual

--------------------------------------------------------------------------------
Adding Geometrical Properties
--------------------------------------------------------------------------------

The network only contain topological and spatial information, so it is necessary to add geometrical information by creating a Geometry object:

.. code-block:: python

    >>> geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)

In this case the `StickAndBall` class was used, which has preset pore-scale models that calculate properties such as diameters and volumes, based on the assumption that the pores are spherical and the throats are cylinders.

--------------------------------------------------------------------------------
Creating Phases
--------------------------------------------------------------------------------

Phases must created to calculate the thermophysical properties of the fluids (and solids) used in the simulations.

    >>> Hg = op.phases.Mercury(network=pn)

OpenPNM includes a few common phases, including  :ref:`air_api` and :ref:`water_api`, but also a wealth of pore-scale models for calculating properties of different phases.


    >>> phys = op.physics.Standard(network=pn, phase=Hg, geometry=geo)

The network can be visualized in `Paraview <http://www.paraview.org>`_ giving the following:

.. image:: http://i.imgur.com/GbUNy0b.png
   :width: 500 px
   :align: center

The drainage curve can be visualized with ``MIP.plot_drainage_curve()`` giving something like this:

.. image:: http://i.imgur.com/ZxuCict.png
   :width: 500 px
   :align: center

A collection of examples has been started as a new Github repository: `OpenPNM-Examples <https://www.github.com/PMEAL/OpenPNM-Examples>`_.
