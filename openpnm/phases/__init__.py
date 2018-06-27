r"""
===============================================================================
openpnm.phases
===============================================================================

This module contains the GenericPhase class, plus several subclasses which
are preconfigured to have the properties of specific fluids

"""
from .GenericPhase import GenericPhase
from .Air import Air
from .Water import Water
from .Mercury import Mercury
