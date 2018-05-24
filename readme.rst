.. image:: https://badge.fury.io/py/openpnm.svg
   :target: https://pypi.python.org/pypi/openpnm

.. image:: https://travis-ci.org/PMEAL/OpenPNM.svg?branch=master
   :target: https://travis-ci.org/PMEAL/OpenPNM

.. image:: https://codecov.io/gh/PMEAL/OpenPNM/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/PMEAL/OpenPNM

.. image:: https://readthedocs.org/projects/openpnm/badge/?version=latest
   :target: http://openpnm.readthedocs.org/

###############################################################################
Overview of OpenPNM
###############################################################################

*OpenPNM* is an open source project aiming to provide porous media researchers with a ready-made framework for performing a wide range of pore network simulations.  The main features and capabilities of OpenPNM are:

.. list-table::

    * - **Defines a universal means of representing any network topology based on a sparse representation of the adjacency matrix**
        - Includes network generators for creating cubic or random networks with arbitrary connectivity
    * - **Provides a set of tools for querying, inspecting, and manipulating topology**
        - Including finding neighboring pores, labeling specific locations, adding or removing pores and throats, joining networks, subdividing and merging pores to create multiscale models, and much more
    * - **Stores pore and throat property data in vectorized format**
        - Allows for fast calculations even on large networks
        - Supports the familiar and advanced array access features such as direct indexing, slicing, Boolean masking, etc.
    * - **A mechanism for calculating the pore-scale properties that define the geometrical (i.e. pore radius), thermophysical (i.e. viscosity), and physics (i.e. hydraulic conductance) properties of the simulation**
        - The interdependence of some properties on other properties is naturally included so values can be regenerated when changes occur (i.e. viscosity can be updated when temperature changed)
        - This mechanism was designed to allow users to easily create new customized pore-scale models suitable for their specific domain
        - A wide assortment of pore-scale transport parameter, pore size calculations, and thermophysical property models are included
    * - **A suite of algorithms for performing network simulations**
        - Including invasion percolation, capillary drainage, mass diffusion, permeability and so on.
        - All transport algorithms can be run transiently and with multiple reactions
    * - **Supports saving, loading, importing and exporting data in numerous formats**
        -  Allows importing networks generated or extracted by other codes, as well as exporting data for post-processing and visualization
        - Saving and loading of simulations allows for batch processing of simulations to be analyzed at a later point

===============================================================================
Installation
===============================================================================

OpenPNM can be installed from the Python Package index using:

.. code-block::

   pip install openpnm

Or the source code can be downloaded from `Github <https://github.com/pmeal/OpenPNM/>`_ and installed by running:

.. code-block::

   python setup.py

===============================================================================
Example Usage
===============================================================================

The following code block illustrates how to use OpenPNM to perform a mercury intrusion porosimetry simulation in just 10 lines:

.. code-block:: python

    >>> import OpenPNM as op
    >>> pn = op.Network.Cubic(shape=[10, 10, 10], spacing=0.0001)
    >>> geo = op.Geometry.Stick_and_Ball(network=pn, pores=pn.Ps,
    ...                                  throats=pn.Ts)
    >>> Hg = op.Phases.Mercury(network=pn)
    >>> Air = op.Phases.Air(network=pn)
    >>> phys = op.Physics.Standard(network=pn, phase=Hg, pores=pn.Ps,
    ...                            throats=pn.Ts)
    >>> MIP = op.Algorithms.Drainage(network=pn)
    >>> MIP.setup(invading_phase=Hg, defending_phase=Air)
    >>> MIP.set_inlets(pores=pn.pores(['top', 'bottom']))
    >>> MIP.run()

The network can be visualized in `Paraview <http://www.paraview.org>`_ giving the following:

.. image:: https://i.imgur.com/mSDrIBOm.png

The drainage curve can be visualized with ``MIP.plot_drainage_curve()`` giving something like this:

.. image:: https://i.imgur.com/1C2uXt9m.png

A collection of examples is available as a separate Github repository: `OpenPNM-Examples <https://www.github.com/PMEAL/OpenPNM-Examples>`_.

===============================================================================
Release Management and Versioning
===============================================================================

OpenPNM uses `Semantic Versioning <http://semver.org>`_ (i.e. X.Y.Z) to label releases.  All major and minor versions (X.Y.z) are available on `PyPI <https://pypi.python.org/pypi>`_, but bugfixes (x.y.Z) updates are not generally pushed unless the bug is particularly egregious.

OpenPNM uses the `Github Flow <http://scottchacon.com/2011/08/31/github-flow.html>`_ system of Git branching. Any code added to master is done via Pull Requests (PRs).  When new PRs are merged into the master branch, they are generally *not* given a new version number. Once enough new features have been added, or a sufficient amount of time has passed, the minor release number (x.Y.z) will be incremented. Any code added to the Master branch *between* incrementing the version number is subject to change, but once a version number has been tagged to code can be considered up-to-date and stable.

OpenPNM depends on several other packages widely known as the `Scipy Stack <https://www.scipy.org/stackspec.html>`_.  It is our policy to always support the latest version of all these packages and their dependencies.

The main developer for this project is Prof. Jeff Gostick (jgostick@gmail.com).

===============================================================================
Licence and Citation
===============================================================================

OpenPNM is free to use and is offered under the permissive `MIT License <http://opensource.org/licenses/MIT>`_.

If you do use OpenPNM in an academic work, the developers ask that you cite the following paper, which outlines the design principles and general uses of OpenPNM:

::

    Gostick et al. OpenPNM: A pore network modeling package. Computing in Science & Engineering. 18(4), p60-74.

A link to this article can be found `here <http://doi.org/10.1109/MCSE.2016.49>`_.
