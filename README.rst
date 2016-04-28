.. image:: https://badge.fury.io/py/openpnm.svg
   :target: https://pypi.python.org/pypi/openpnm

.. image:: https://travis-ci.org/PMEAL/OpenPNM.svg?branch=master
   :target: https://travis-ci.org/PMEAL/OpenPNM

.. image:: https://codecov.io/gh/PMEAL/OpenPNM/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/PMEAL/OpenPNM

.. image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Join the chat at https://gitter.im/PMEAL/OpenPNM
   :target: https://gitter.im/PMEAL/OpenPNM?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge

###############################################################################
Overview of OpenPNM
###############################################################################

*OpenPNM* is an open source project aiming to provide porous media researchers with a ready-made framework for performing a wide range of pore network simulations.  The main features and capabilities of OpenPNM are:

    * Defines a universal means of representing any network topology based on a sparse representation of the adjacency matrix

    * Provides a set of tools for querying, inspecting, and manipulating topology

      * Including finding neighboring pores, labeling specific locations, adding or removing pores and throats, joining networks, subdividing and merging pores to create multiscale models, and much more

    * Includes network generators for creating cubic or random networks with arbitrary connectivity

    * Stores pore and throat property data in vectorized format

      * Allows for fast calculations even on large networks

      * Supports the familiar and advanced array access features such as direct indexing, slicing, Boolean masking, etc.

    * A sophisticated mechanism for calculating the pore-scale properties that define the geometrical (i.e. pore radius), thermophysical (i.e. viscosity), and physics (i.e. hydraulic conductance) properties of the simulation

      * The interdependence of some properties on other properties is naturally included so values can be regenerated when changes occur (i.e. viscosity can be updated when temperature changed)

      * This mechanism was designed to allow users to easily create new customized pore-scale models suitable for their specific domain

      * A wide assortment of pore-scale transport parameter, pore size calculations, and thermophysical property models are included

    * A suite of algorithms for performing network simulations such as invasion percolation, capillary drainage, mass diffusion, permeability and so on.

    * Supports saving, loading, importing and exporting data in numerous formats

      *  Allows importing networks generated or extracted by other code, as well as exporting data for post-processing and visualization

      * Saving and loading of simulations allows for batch processing of simulations to be analyzed at a later point

===============================================================================
Installation
===============================================================================

OpenPNM can be install from the Python Package index using:

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

.. image:: http://i.imgur.com/GbUNy0bm.png

The drainage curve can be visualized with ``MIP.plot_drainage_curve()`` giving something like this:

.. image:: http://i.imgur.com/ZxuCictm.png

A collection of examples is available as a separate Github repository: `OpenPNM-Examples <https://www.github.com/PMEAL/OpenPNM-Examples>`_.

===============================================================================
Related Links
===============================================================================

.. list-table::

    * - OpenPNM Homepage
        - http://openpnm.org
    * - Github is used to host the code
        - https://www.github.com/PMEAL/OpenPNM
    * - Github is also used as the project's issue and bug tracker
        - https://www.github.com/PMEAL/OpenPNM/issues
    * - A collection of examples using OpenPNM is available in a separate repository
        - https://www.github.com/PMEAL/OpenPNM-Examples
    * - Gitter is used to help users with questions about using the code
        - https://gitter.im/PMEAL/OpenPNM
    * - Scipy is a major component of OpenPNM
        - http://www.scipy.org
    * - Anaconda is the most general way to setup a numerical Python environment
        - https://www.continuum.io/downloads
    * - WinPython is a slightly easier way to use numerical Python on Windows
        - https://github.com/winpython/winpython
    * - Spyder is the recommended IDE when working with OpenPNM
        - https://github.com/spyder-ide/spyder
    * - Paraview is suggested for visualizing OpenPNM data
        - http://www.paraview.org
    * - OpenPNM is offered under an MIT License
        - http://opensource.org/licenses/MIT
