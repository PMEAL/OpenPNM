|
.. _front_page:


.. image:: https://travis-ci.org/PMEAL/OpenPNM.svg?branch=master
   :target: https://travis-ci.org/PMEAL/OpenPNM

.. image:: https://img.shields.io/codecov/c/github/PMEAL/OpenPNM.svg?style=flat
   :target: https://codecov.io/github/PMEAL/OpenPNM


###############################################################################
Overview of OpenPNM
###############################################################################

*OpenPNM* is an open source project aiming to provide porous media researchers with a ready-made framework for performing a wide range of pore network simulations.

===============================================================================
Documentation
===============================================================================

.. toctree::
    :maxdepth: 3

    userguide/index.rst
    docOpenPNM.rst

===============================================================================
Example Usage
===============================================================================

The following code block illustrates how to use OpenPNM to perform a mercury intrusion porosimetry simulation in just 10 lines:

.. code-block:: python

    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[10, 10, 10], spacing=0.0001)
    >>> geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
    >>> Hg = op.phases.Mercury(network=pn)
    >>> Air = op.phases.Air(network=pn)
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
