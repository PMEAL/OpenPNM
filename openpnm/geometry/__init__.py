r"""

**openpnm.geometry**

----

The ``geometry`` module contains the ``GenericGeometry`` class, and an
assortment of subclasses that implement specific pore-scale geometrical models.


"""
from . import collections
from .GenericGeometry import GenericGeometry
from .Imported import Imported
from .SpheresAndCylinders import SpheresAndCylinders
from .CirclesAndRectangles import CirclesAndRectangles
from .ConesAndCylinders import ConesAndCylinders
from .PyramidsAndCuboids import PyramidsAndCuboids
from .TrapezoidsAndRectangles import TrapezoidsAndRectangles
from .CubesAndCuboids import CubesAndCuboids
from .SquaresAndRectangles import SquaresAndRectangles
