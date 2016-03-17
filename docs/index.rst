###############################################################################
OpenPNM: Pore Network Modeling in Python
###############################################################################

*OpenPNM* is an open source project aiming to provide the scientific and engineering community with a ready-made framework for performing pore network simulations.  Pore network models have been used for decades among porous media researchers as a way to efficiently model multiphase transport and capillary phenomena in porous media.  The *OpenPNM* framework supplies a means of storing network data, representing network connections, performing queries on the network, adding/removing pores and so on.  OpenPNM also provides a suite of algorithms for performing simulations on the network such as invasion percolation, capillary drainage, mass diffusion, permeability and so on.

The following code snippet illustrates how to use OpenPNM to perform a mercury intrusion porosimetry simulation in just 10 lines:

.. code-block:: python

    >>> import OpenPNM as op
    >>> pn = op.Network.Cubic(shape=[10, 10, 10], spacing=0.0001)
    >>> geo = op.Geometry.Stick_and_Ball(network=pn, pores=pn.Ps, throats=pn.Ts)
    >>> # Key size distributions can be visualized with geo.plot_histograms()
    >>> Hg = op.Phases.Mercury(network=pn)
    >>> Air = op.Phases.Air(network=pn)
    >>> phys = op.Physics.Standard(network=pn, phase=Hg, pores=pn.Ps, throats=pn.Ts)
    >>> MIP = op.Algorithms.Drainage(network=pn)
    >>> MIP.setup(invading_phase=Hg, defending_phase=Air)
    >>> MIP.set_inlets(pores=pn.pores(['top', 'bottom']))
    >>> MIP.run()
    >>> # The drainage curve can be visualized with MIP.plot_drainage_curve()



+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Links
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  .. list-table::

      * - OpenPNM homepage
	      - http://www.openpnm.
      * - Github is used to host the code
          - http://www.github.com/PMEAL/OpenPNM
      * - Github is also used as the project's issue and bug tracker
          - http://www.github.com/PMEAL/OpenPNM/issues
      * - Spyder is the recommended IDE when working with OpenPNM
          - https://github.com/spyder-ide/spyder
      * - Paraview is suggested for visualizing OpenPNM data
          - http://paraview.org
      * - Scipy is a major component of OpenPNM
          - http://www.scipy
      * - Anaconda is the recommended way to setup a numerical Python environment
          - https://www.continuum.io/downloads
      * - OpenPNM is offered under an MIT License
          - http://opensource.org/licenses/MIT

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Documentation
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. toctree::
   :maxdepth: 2

   userguide/index.rst
   devguide/index.rst
   docOpenPNM.rst
   docExtraSidebar.rst
