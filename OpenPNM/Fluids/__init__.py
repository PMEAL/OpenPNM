r"""
*******************************************************************************
:mod:`OpenPNM.Fluids`: Fluid Property Estimation Methods
*******************************************************************************

.. module:: OpenPNM.Fluids

Contents
--------
This submodule contains methods for estimating fluid properties

.. note::
    none

.. autoclass:: GenericFluid
   :members:
   :undoc-members:
   :private-members:
   :show-inheritance:
   
.. autoclass:: Water
   :members:
   :undoc-members:
   :private-members:
   :show-inheritance:
   
.. autoclass:: Air
   :members:
   :undoc-members:
   :private-members:
   :show-inheritance:
   
.. autoclass:: Viscosity
   :members:
   :undoc-members:
   :private-members:
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