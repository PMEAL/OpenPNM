r"""
*********************************************************************************
:mod:`OpenPNM.Geometry` -- All classes related the creation of network geometry
*********************************************************************************

.. module:: OpenPNM.Geometry

Contents
--------
The OpenPNM package imports all the functions from the top level modules. 

.. note::
    n/a
 
Import
------
>>> import OpenPNM as PNM
>>> tmp=PNM.Geometry.GenericGeometry()


Submodules
----------
::

 None                            --- No subpackages at the moment
 
Classes
-------
    
.. autoclass:: GenericGeometry
   :members:
   :undoc-members:
   :show-inheritance:
       
.. autoclass:: Cubic
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: Delaunay
   :members:
   :undoc-members:
   :show-inheritance:
   
.. autoclass:: Template
   :members:
   :undoc-members:
   :show-inheritance:
   
"""

import scipy as sp
import numpy as np

from __GenericGeometry__ import GenericGeometry
from __Cubic__ import Cubic
from __Delaunay__ import Delaunay
from __Template__ import Template
from __Import__ import MatFile
