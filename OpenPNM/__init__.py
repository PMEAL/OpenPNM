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
   * - :mod:`OpenPNM.BAS`
     - common utilities and classes used by most of the of the modules
   * - :mod:`OpenPNM.NET`
     - Storage and manipulations of network topoologies and data stored on them.
   * - :mod:`OpenPNM.GEN`
     - Generators for pore networks. (Random cubic, image based, Voronoi). Should also contain
       a mapper of the pore network back on the segmented image.
   * - :mod:`OpenPNM.ALG`
     - Module containing all algorithmic classes for networks.
   * - `IO`
     - Input output routines
   * - `VISU`
     - Mayavi-based post-processing modules (`postproc.py`)   
   * - `interactive/`
     - setup of IPython-based shell `isfepy`
   * - `linalg/`
     - linear algebra functions not covered by NumPy and SciPy


 
 
Utility tools
-------------
::

 TODO                --- Todo
 
 
Import
------
>>> import OpenPNM as PNM

Inheritance Diagram
--------------------


.. inheritance-diagram:: OpenPNM.NET.GenericNetwork

Package Documentation
---------------------

.. automodule:: BAS  
   :members:
   :undoc-members:
   :show-inheritance:

      
.. automodule:: NET
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: GEN
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: ALG
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

# __all__ = ['BAS']



import BAS
import NET
import GEN
import IO
import PHYS
import ALG
import IMG
import GUI
import VIS
#
# from PNMbase import *
# from PNMalgorithms import *
# from PNMgenerators import *
# from PNMIO import *
# from PNMparameters import *
# from PNMtransport import *



