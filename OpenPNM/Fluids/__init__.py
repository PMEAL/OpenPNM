r"""
*******************************************************************************
:mod:`OpenPNM.Fluids` -- Fluid Property Estimation Methods
*******************************************************************************

.. module:: OpenPNM.Fluids

Contents
--------
This module contains methods for estimating fluid properties

.. autoclass:: GenericFluid
   :members:
   :undoc-members:
   :show-inheritance:
   
.. autoclass:: Water
   :members:
   :undoc-members:
   :show-inheritance:
   
.. autoclass:: Air
   :members:
   :undoc-members:
   :show-inheritance:
   
.. automodule:: OpenPNM.Fluids.contact_angle
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Fluids.diffusivity
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Fluids.electrical_conductivity
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Fluids.molar_density
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Fluids.molar_mass
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Fluids.surface_tension
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Fluids.thermal_conductivity
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Fluids.vapor_pressure
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Fluids.viscosity
   :members:
   :undoc-members:
   :show-inheritance:

"""

from .__GenericFluid__ import GenericFluid
from .__Water__ import Water
from .__Air__ import Air
from . import viscosity
from . import molar_density
from . import molar_mass
from . import diffusivity
from . import surface_tension
from . import vapor_pressure
from . import contact_angle
from . import electrical_conductivity
from . import thermal_conductivity