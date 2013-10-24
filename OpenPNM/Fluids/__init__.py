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

from __GenericFluid__ import GenericFluid
from __Water__ import Water
from __Air__ import Air
import Viscosity
import MolarDensity
import Diffusivity
import SurfaceTension
import VaporPressure