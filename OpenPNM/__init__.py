# -*- coding: utf-8 -*-
# Author: Andreas Putz
# Copyright (c) 2013, OpenPNM
# License: TBD.

r"""
##################################################################################
:mod:`OpenPNM` --  A scientific pore network calculator for porous transport media
##################################################################################
.. module:: OpenPNM
    :platform: Linux, Windows

Documentation is available in the docstrings and in ths sphinx documentation.

Contents
--------
The OpenPNM package imports all the functions from the top level modules.


Subpackages
-----------


.. list-table:: OpenPNM submodule structure.
   :widths: 10 80 
   :header-rows: 1

   * - Name
     - Description
   * - :mod:`OpenPNM.Utilities`
     - common utilities and classes used by most of the of the modules
   * - :mod:`OpenPNM.Network`
     - Storage and manipulations of network topoologies and data stored on them.
   * - :mod:`OpenPNM.Geometry`
     - Geometry for pore networks. (Random cubic, image based, Voronoi). Should also contain
       a mapper of the pore network back on the segmented image.
   * - :mod:`OpenPNM.Algorithsm`
     - Module containing all algorithmic classes for networks.
   * - `Visualization`
     - vtk-based post-processing modules (`postproc.py`)   


 
 
Utility tools
-------------
::

 TODO                --- Todo
 
 
Import
------
>>> import OpenPNM as PNM

Inheritance Diagram
--------------------


.. inheritance-diagram:: OpenPNM.Network.GenericNetwork

Package Documentation
---------------------

.. automodule:: Base
   :members:
   :undoc-members:
   :show-inheritance:

      
.. automodule:: Network
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Geometry
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Algorithms
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Physics
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: GUI
   :members:
   :undoc-members:
   :show-inheritance:
"""

__version__ = '0.0.1'

__requires__ = [
    'scipy',
    'numpy',
]

__extras_require__ = {
    'app': [
        'envisage',
    ],
}

# __all__ = ['Base']


import Utilities
import Network
import Geometry
import Algorithms
import Visualization
import Physics
import GUI





