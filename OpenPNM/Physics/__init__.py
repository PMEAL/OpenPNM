r"""
===============================================================================
:mod:`OpenPNM.Physics` -- Pore Scale Physics Models
===============================================================================

.. module:: OpenPNM.Physics

Contents
--------
This submodule contains all pore scale physics models applied to a pore network.


.. autoclass:: GenericPhysics
   :members:
   :undoc-members:
   :show-inheritance:
   
.. automodule:: OpenPNM.Physics.capillary_pressure
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Physics.diffusive_conductance
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Physics.electronic_conductance
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Physics.hydraulic_conductance
   :members:
   :undoc-members:
   :show-inheritance:
   
.. automodule:: OpenPNM.Physics.thermal_conductance
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
