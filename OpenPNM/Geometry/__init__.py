r"""
###############################################################################
:mod:`OpenPNM.Geometry`: Classes related to the creation of pore-scale geometry
###############################################################################

Contents
--------
GenericGeometry: The init statement performs the important tasks of registering
itself with the network and creating data dictionaries of the necessary size

Subclasses: These can be custom made by users to represent specfic geometries.
It is necessary for their init's to call the init of the GenericGeometry class
in order to be properly instantiated.  OpenPNM comes with a few basic
pre-written materials.  New classes added to this directory will be
automatically imported and available.

Classes
-------

.. autoclass:: GenericGeometry
   :members:

.. autoclass:: Stick_and_Ball
   :members:

.. autoclass:: Cube_and_Cuboid
   :members:

.. autoclass:: Boundary
   :members:

.. autoclass:: Voronoi
   :members:

"""

from .__GenericGeometry__ import GenericGeometry
from .__Cube_and_Cuboid__ import Cube_and_Cuboid
from .__Stick_and_Ball__ import Stick_and_Ball
from .__TestGeometry__ import TestGeometry
from .__Boundary__ import Boundary
from .__Voronoi__ import Voronoi
from .__Toray090__ import Toray090
from .__SGL10__ import SGL10
from . import models
