|
.. _front_page:

.. image:: https://badge.fury.io/py/openpnm.svg
   :target: http://badge.fury.io/py/openpnm

.. image:: https://travis-ci.org/PMEAL/OpenPNM.svg?branch=master
   :target: https://travis-ci.org/PMEAL/OpenPNM

.. image:: https://img.shields.io/codecov/c/github/PMEAL/OpenPNM.svg?style=flat
   :target: https://codecov.io/github/PMEAL/OpenPNM

.. image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Join the chat at https://gitter.im/PMEAL/OpenPNM
   :target: https://gitter.im/PMEAL/OpenPNM?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge

###############################################################################
Overview of OpenPNM
###############################################################################

*OpenPNM* is an open source project aiming to provide porous media researchers with a ready-made framework for performing a wide range of pore network simulations.

.. list-table:: **Summary of key capabilities offered by OpenPNM**

    * - **Defines a universal means of representing any network topology**
      - Based on a sparse representation of the adjacency matrix using principles from graph theory.
    * - **Provides a set of tools for querying, inspecting, and manipulating topology**
      - Including finding neighboring pores, labeling specific locations, adding or removing pores and throats, joining networks, subdividing and merging pores to create multiscale models, and much more.
    * - **Able to generate various network topologies**
      - Includes network generators for creating cubic or random networks with arbitrary connectivity.
    * - **Stores pore and throat property data in vectorized format**
      - Allows for fast calculations even on large networks, and support for the familiar and advanced array access features such as direct indexing, slicing, Boolean masking, etc.
    * - **Includes a sophisticated mechanism for calculating the pore-scale properties**
      - A wide assortment of pore-scale transport parameters, pore size calculations, and thermophysical property models are included, and new models can easily be created by users for their specific problem.
    * - **Ships with a suite of algorithms for performing network simulations**
      - Includes invasion percolation, capillary drainage, mass diffusion, permeability and so on.
    * - **Supports saving, loading, importing and exporting data in numerous formats**
      - Allows importing networks generated or extracted by other code, and exporting data for post-processing and visualization, as well as a native format for saving and loading complete simulations for future analysis.

===============================================================================
Documentation
===============================================================================

.. toctree::
    :maxdepth: 3

    userguide/index.rst
    examples.rst
    docOpenPNM.rst

===============================================================================
Example Usage
===============================================================================

The following code block illustrates how to use OpenPNM to perform a mercury intrusion porosimetry simulation in just 10 lines:

.. code-block:: python

    >>> import OpenPNM as op
    >>> pn = op.Network.Cubic(shape=[10, 10, 10], spacing=0.0001)
    >>> geo = op.Geometry.Stick_and_Ball(network=pn, pores=pn.Ps, throats=pn.Ts)
    >>> Hg = op.Phases.Mercury(network=pn)
    >>> Air = op.Phases.Air(network=pn)
    >>> phys = op.Physics.Standard(network=pn, phase=Hg, pores=pn.Ps, throats=pn.Ts)
    >>> MIP = op.Algorithms.Drainage(network=pn)
    >>> MIP.setup(invading_phase=Hg, defending_phase=Air)
    >>> MIP.set_inlets(pores=pn.pores(['top', 'bottom']))
    >>> MIP.run()

The network can be visualized in `Paraview <http://www.paraview.org>`_ giving the following:

.. image:: http://i.imgur.com/GbUNy0b.png
   :width: 500 px
   :align: center

The drainage curve can be visualized with ``MIP.plot_drainage_curve()`` giving something like this:

.. image:: http://i.imgur.com/ZxuCict.png
   :width: 500 px
   :align: center

A collection of examples has been started as a new Github repository: `OpenPNM-Examples <https://www.github.com/PMEAL/OpenPNM-Examples>`_.
