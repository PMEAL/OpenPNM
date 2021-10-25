r"""

**openpnm.geometry**

----

The ``geometry`` module contains the ``GenericGeometry`` class, and an
assortment of subclasses that implement specific pore-scale geometrical models.

----

**The GenericGeometry Class**

Geometry objects (as well as Physics objects) are ``Subdomain`` subclasses,
which allow them to be assigned to subset of the full domain (although this is
not alway necessary).  This functionality was added so that networks with
distinct regions could be modelled by giving each region its own Geometry
with unique models (e.g. to give a bi-modal pore size distribution).

----

**Library of Preconfigured Geometry Classes**

This module contains a small selection of Geometry classes that are
pre-configured with a selection of pore-scale models.  These classes provide
a good starting point, but generally the choice of models and parameters used
will be specific to each problem and must be designed by the user.

The ``StickAndBall`` class, as it's name suggests assumes spherical pores and
cylindrical throats.  Pore sizes are assigned by finding the largest sphere
that can fit at each site (this will be dictated by the lattice spacing used
when generating the Network), then scaling that value by a random number
between 0 and 0.1.  Throat diameters are taken as half the size of the smaller
of it's two neighbors.  All other properties are calculated using the geometry
of spheres and throats.

The table belows shows the specific models used on ``StickAndBall``:

+----+----------------------+------------------+--------------------------+
| #  | Property Name        | Parameter        | Value                    |
+====+======================+==================+==========================+
| 1  | pore.seed            | model:           | random                   |
+----+----------------------+------------------+--------------------------+
| 2  | pore.max_size        | model:           | largest_sphere           |
+----+----------------------+------------------+--------------------------+
| 3  | pore.diameter        | model:           | product                  |
+----+----------------------+------------------+--------------------------+
| 4  | pore.area            | model:           | sphere                   |
+----+----------------------+------------------+--------------------------+
| 5  | pore.volume          | model:           | sphere                   |
+----+----------------------+------------------+--------------------------+
| 6  | throat.max_size      | model:           | from_neighbor_pores      |
+----+----------------------+------------------+--------------------------+
| 7  | throat.diameter      | model:           | scaled                   |
+----+----------------------+------------------+--------------------------+
| 8  | throat.length        | model:           | piecewise                |
+----+----------------------+------------------+--------------------------+
| 9  | throat.surface_area  | model:           | cylinder                 |
+----+----------------------+------------------+--------------------------+
| 10 | throat.volume        | model:           | cylinder                 |
+----+----------------------+------------------+--------------------------+
| 11 | throat.area          | model:           | cylinder                 |
+----+----------------------+------------------+--------------------------+

----

**Customizing a Preconfigured Geometry Instance**

Perhaps the ``StickAndBall`` class is almost suitable but you wish to decrease
the pores sizes.  The following example illustrates how to alter the
``'pore.size'`` model accordingly:

>>> import openpnm as op
>>> pn = op.network.Cubic([5, 5, 5])
>>> geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)

We can reach into the ``models`` attribute and change the parameters of any
model as follows:

>>> max(geo['pore.diameter']) < 1.0  # Confirm largest pore is less than 1.0
True
>>> geo.models['pore.seed']['num_range'] = [0.2, 0.7]
>>> geo.regenerate_models()  # Must regenerate all models
>>> max(geo['pore.diameter']) < 0.7  # Largest pore is now less than 0.7
True

This example illustrated that you can change one property ('pore.seed') and
that change can be cascaded to all dependent properties ('pore.diameter').

"""
from .GenericGeometry import GenericGeometry
from .StickAndBall import StickAndBall
from .StickAndBall2D import StickAndBall2D
from .Boundary import Boundary
from .Boundary2D import Boundary2D
from .Imported import Imported
from .SpheresAndCylinders import SpheresAndCylinders
from .CirclesAndRectangles import CirclesAndRectangles
from .ConesAndCylinders import ConesAndCylinders
from .PyramidsAndCuboids import PyramidsAndCuboids
from .TrapezoidsAndRectangles import TrapezoidsAndRectangles
from .CubesAndCuboids import CubesAndCuboids
from .SquaresAndRectangles import SquaresAndRectangles
