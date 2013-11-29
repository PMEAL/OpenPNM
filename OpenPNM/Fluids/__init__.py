r"""
*******************************************************************************
:mod:`OpenPNM.Fluids`: FluidProperties
*******************************************************************************

.. module:: OpenPNM.Fluids

Contents
--------
This submodule contains methods for estimating fluid properties

.. note::
    none


Submodules
----------
::

 None                            --- No subpackages at the moment


"""

from .__GenericFluid__ import GenericFluid
from .__Water__ import Water
from .__Air__ import Air
from . import Viscosity
from . import MolarDensity
from . import Diffusivity
from . import SurfaceTension
from . import VaporPressure
from . import ContactAngle
from . import ElectricalConductivity
from . import ThermalConductivity