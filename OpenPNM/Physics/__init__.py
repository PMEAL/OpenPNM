r"""
**************************************************************************
:mod:`OpenPNM.PHYS`: Pore Scale Physics for Networks
**************************************************************************

.. module:: OpenPNM.Physics

Contents
--------
This submodule contains all pore scale physics models applied to a pore network.

.. note::
    The algorithms take a basenet as an argument in the constructor, this
    seems to initialize a full object. Causing a lot of garbage to be written.
 
Import
------
>>> import OpenPNM as PNM
>>> tmp=PNM.Physics.GenericPhysics()


Submodules
----------
::

 None                            --- No subpackages at the moment

.. autoclass:: Physics.GenericPhysics
   :members:
   :undoc-members:
   :show-inheritance:

    
"""

from __GenericPhysics__     import GenericPhysics
from __CapillaryPressure__  import CapillaryPressure
from __MassTransport__      import MassTransport
from __HeatConduction__     import HeatConduction
from __FluidFlow__          import FluidFlow
from __Multiphase__         import Multiphase
from __Thermo__             import Thermo