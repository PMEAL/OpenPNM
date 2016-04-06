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
OpenPNM Documentation
###############################################################################

.. toctree::
    :maxdepth: 3

    userguide/index.rst

===============================================================================
Overview of OpenPNM
===============================================================================

*OpenPNM* is an open source project aiming to provide the scientific and engineering community with a ready-made framework for performing pore network simulations.  Pore network models have been used for decades among porous media researchers as a way to efficiently model multiphase transport and capillary phenomena in porous media.  The *OpenPNM* framework supplies a means of representing network connections, storing geometrical data, c.

OpenPNM also provides a suite of algorithms for performing simulations on the network such as invasion percolation, capillary drainage, mass diffusion, permeability and so on.

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
    * - Spyder is the recommended IDE when working with OpenPNM
        - https://github.com/spyder-ide/spyder
    * - Paraview is suggested for visualizing OpenPNM data
        - http://www.paraview.org
    * - Scipy is a major component of OpenPNM
        - http://www.scipy.org
    * - Anaconda is the most general way to setup a numerical Python environment
        - https://www.continuum.io/downloads
    * - WinPython is a slightly easier way to use numerical Python on Windows
        - https://github.com/winpython/winpython
    * - OpenPNM is offered under an MIT License
        - http://opensource.org/licenses/MIT
