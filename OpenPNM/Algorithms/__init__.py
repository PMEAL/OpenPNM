# Author: Andreas Putz
# Copyright (c) 2013, OpenPNM
# License: TBD.
r"""
**************************************************************************
:mod:`OpenPNM.Algorithms`: Algorithms on Networks
**************************************************************************

.. module:: OpenPNM.Algorithms

Contents
--------
This submodule contains all algorithms actiong on a pore network.

.. note::
    The algorithms take a basenet as an argument in the constructor, this
    seems to initialize a full object. Causing a lot of garbage to be written.
 
Import
------
>>> import OpenPNM as PNM
>>> tmp=PNM.Algorithms.GenericAlgorithm()


Submodules
----------
::

 None                            --- No subpackages at the moment

.. autoclass:: Algorithms.GenericAlgorithm
   :members:
   :undoc-members:
   :show-inheritance:
       
.. autoclass:: Algorithms.InvasionPercolation
   :members:
   :undoc-members:
   :show-inheritance:
   
.. autoclass:: Algorithms.OrdinaryPercolation
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: Algorithms.FickianDiffusion
   :members:
   :undoc-members:
   :show-inheritance:
    
"""

from __GenericAlgorithm__ import GenericAlgorithm
from __InvasionPercolation__ import InvasionPercolation
from __OrdinaryPercolation__ import OrdinaryPercolation
from __FickianDiffusion__ import FickianDiffusion 