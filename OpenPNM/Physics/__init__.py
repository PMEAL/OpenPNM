r"""
*******************************************************************************
:mod:`OpenPNM.Physics`: Pore Scale Physics Models
*******************************************************************************

.. module:: OpenPNM.Physics

Contents
--------
This submodule contains all pore scale physics models applied to a pore network.

.. note::
    none


Submodules
----------
::

 None                            --- No subpackages at the moment

.. automodule:: Physics.capillary_pressure
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Physics.hydraulic_conductance
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Physics.diffusive_conductance
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Physics.electronic_conductance
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Physics.thermal_conductance
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Physics.multi_phase
   :members:
   :undoc-members:
   :show-inheritance:

"""
from .__GenericPhysics__ import GenericPhysics
from . import electronic_conductance
from . import capillary_pressure
from . import diffusive_conductance
from . import thermal_conductance
from . import hydraulic_conductance
from . import multi_phase