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


"""

from .__GenericPhysics__ import GenericPhysics
from . import electronic_conductance
from . import capillary_pressure
from . import diffusive_conductance
from . import thermal_conductance
from . import hydraulic_conductance
