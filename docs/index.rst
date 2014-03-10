.. OpenPNM documentation master file, created by
   sphinx-quickstart on Thu Nov 29 11:31:33 2012.
   
*******************************************************************************
OpenPNM: Pore Network Modeling in Python
*******************************************************************************

*OpenPNM* is an open source pore network modeling package initiated by researchers at McGill University and the University of Toronto and sponsored by the Automotive Fuel Cell Cooperation.  The aim of this package is to provide the scientific and engineering community with a ready-made framework for performing pore network simulations.  Pore network models have been used for decades among porous media researchers as a way to efficiently and accurately model multiphase transport and capillary phenomena.  Until now all researchers in this field must build their code from scratch since there is no mainstream commercial offering, as there is for CFD based modeling (ie. COMSOL and FLUENT).  OpenPNM aims to remedy this situation by providing a fast, efficient framework for working with pore network models of arbitrary size and dimensionality.

*OpenPNM* is coded in Python relying heavily on `Scipy <http://www.scipy.org>`_.  Python is a free and open source programming language that compares very closely with Matlab.  The *OpenPNM* framework supplies a means of storing network data, representing network connections, performing queries on the network, adding/removing pores and so on.  OpenPNM also provides a suite of algorithms for performing simulations on the network such as invasion percolation, capillary drainage, mass diffusion, permeability and so on.  Each algorithm is stored as a separate object so new algorithms can be contributed by users as standalone objects without having to alter the basic network.

Working with *OpenPNM* is most effectively accomplished using scripts, but a GUI is also being developed to provide immediate access to the basic features of the software.  At present Paraview is used for visualization purposes.  *OpenPNM* outputs a VTK file  which can be imported into Paraview as a *polydata* file format.  Paraview is perfectly suited to visualizing the results of *OpenPNM*.

*OpenPNM* is not yet ready for distribution.  It is envisioned that version 1.0 will be available around before Summer 2014.  The source code is hosted on Github for those wishing to take a closer look.


**License:** `MIT <http://opensource.org/licenses/MIT>`_

Links
=====

  .. list-table::

      * - OpenPNM homepage
	- http://www.openpnm.org
      * - GIT project site
        - https://github.com/PMEAL/OpenPNM

.. _documentation:

Documentation
=============

Contents:

.. toctree::
   :maxdepth: 2

   userguide/index.rst
   devguide/index.rst
   docOpenPNM.rst


Index and Tables
================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

