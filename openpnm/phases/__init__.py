r"""

**openpnm.phases**

----

This module contains the GenericPhase class, plus several subclasses which
are preconfigured to have the properties of specific fluids

----

**The Phase Library**

OpenPNM include a few Phase subclasses that contain a suite of pre-configured
models that predict the thermophysical properties of certain common phases.
Specifically, these are: ``Air``, ``Water``, and ``Mercury``.

**The GenericPhase**

Phase objects apply to the entire domain.

"""
from .GenericPhase import GenericPhase
from .Air import Air
from .Water import Water
from .Mercury import Mercury
