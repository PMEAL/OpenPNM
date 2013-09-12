.. FCSTpython documentation master file, created by
   sphinx-quickstart on Thu Nov 29 11:31:33 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
**************************************
OpenPNM: Pore Network Models in Python
**************************************




*OpenPNM* is an open source pore network modeling package initiated by researchers at McGill University and the University of Toronto and sponsored by the Automotive Fuel Cell Cooperation.  The aim of this package is to provide the scientific and engineering community with a ready-made framework for performing pore network simulations.  Pore network models have been used for decades among porous media researchers as a way to efficiently and accurately model multiphase transport and capillary phenomena.  Until now all researchers in this field must build their code from scratch since there is no mainstream commercial offering, as there is for CFD based modeling (ie. COMSOL and FLUENT).  OpenPNM aims to remedy this situation by providing a fast, efficient framework for working with pore network models of arbitrary size and dimensionality.

*OpenPNM* is coded in Python relying heavily on `Scipy <http://www.scipy.org>`_.  Python is a free and open source programming language that compares very closely with Matlab.  The *OpenPNM* framework supplies a means of storing network data, representing network connections, performing queries on the network, adding/removing pores and so on.  The network data and the functions that work with the network are both stored together in a network object.  OpenPNM also provides a suite of algorithms for performing simulations on the network such as invasion percolation, capillary drainage, mass diffusion, permeability and so on.  Each algorithm is stored as a separate object and it is envisioned that new algorithms can be contributed by the user community as standalone objects without having to alter the basic network storage class.

Working with *OpenPNM* is most effectively accomplished using scripts, but a GUI is also being developed to provide immediate access to the basic features of the software.  At present Paraview is used for visualization purposes.  *OpenPNM* outputs a VTK file  which can be imported into Paraview as a polydata file format.  Paraview is perfectly suited to visualizing the results of *OpenPNM*.

*OpenPNM* is not yet ready for distribution.  It is envisioned that version 1.0 will be available around November 1, 2013.  The source code is hosted on Github for those wishing to take a closer look.


**License:** `MIT <http://opensource.org/licenses/MIT>`_

Links
=====

  .. list-table::

      * - OpenPNM homepage
	- http://www.openpnm.org
      * - GIT project site
        - https://github.com/PMEAL/OpenPNM

	
Applications
============

OpenPNM is being developped for applications in the following fields

- Two phase flow in porous media
- Diffusion problems in porous media

Features
========

- Flexible pore network generation based on
  - propability distributions
  - image processing
  - cubic meshes
  - fuel cell specific geometries
- Efficient topological storage
- Pore network algorithms
  - Percolation
  - Access limited percolation
  - Diffusion

Support
=======

Work on OpenPNM has been supported by the following projects:

.. _documentation:

Documentation
=============

Contents:

.. toctree::
   :maxdepth: 2

   docInstall.rst
   docOpenPNM.rst

.. todolist::







Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

