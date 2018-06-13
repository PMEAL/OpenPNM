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

*OpenPNM* is an open source project aiming to provide porous media researchers with a ready-made framework for performing a wide range of pore network simulations.

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

    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[10, 10, 10], spacing=0.0001)
    >>> geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps,
    ...                                throats=pn.Ts)
    >>> Hg = op.phases.Mercury(network=pn)
    >>> Air = op.phases.Air(network=pn)
    >>> phys = op.physics.Standard(network=pn, phase=Hg, geometry=geo)


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
