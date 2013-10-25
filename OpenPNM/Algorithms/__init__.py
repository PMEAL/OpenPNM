r"""
*******************************************************************************
:mod:`OpenPNM.Algorithms`: Algorithms on Networks
*******************************************************************************

.. module:: OpenPNM.Algorithms

Contents
--------
This submodule contains algorithms for performing simulations on pore networks

.. note::
    n/a
 
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
from __Permeability__ import Permeability 
from __FourierConduction__ import FourierConduction 
from __ElectronConduction__ import ElectronConduction
from __LinearSolver__ import LinearSolver
   