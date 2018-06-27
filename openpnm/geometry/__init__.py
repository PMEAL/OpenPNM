r"""
================================================================================
openpnm.geometry
================================================================================

The ``geometry`` module contains a GenericGeometry class, and an assortment
of subclasses that implement specific pore-scale geometrical models

Geometry objects (as well as Physics objects) are Subdomain subclasses, which
allow them to be assigned to subset of the full domain.

"""
from .GenericGeometry import GenericGeometry
from .StickAndBall import StickAndBall
