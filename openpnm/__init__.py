r"""

::

     o-o                o--o  o   o o   o
    o   o               |   | |\  | |\ /|
    |   | o-o  o-o o-o  o--o  | \ | | o |
    o   o |  | |-' |  | |     |  \| |   |
     o-o  o-o  o-o o  o o     o   o o   o
          |
          o

**OpenPNM**

OpenPNM is a package for performing pore network simulations of transport in
porous materials.

It consists of the following submodules:

+----------------+------------------------------------------------------------+
| Submodule      | Contents and Description                                   |
+================+============================================================+
| ``core``       | Houses the ``Base`` & ``Subdomain`` classes, & model       |
|                | related mixins                                             |
+----------------+------------------------------------------------------------+
| ``network``    | ``GenericNetwork`` class plus various network generators   |
+----------------+------------------------------------------------------------+
| ``geometry``   | ``GenericGeometry`` class plus some subclasses containing a|
|                | predefined set of pore-scale models                        |
+----------------+------------------------------------------------------------+
| ``phases``     | ``GenericPhase`` class plus some subclasses containing     |
|                | predefined models for common fluids like Water             |
+----------------+------------------------------------------------------------+
| ``physics``    | ``GenericPhysics`` class plus some subclasses containing a |
|                | predefined set of pore-scale models                        |
+----------------+------------------------------------------------------------+
| ``algorithms`` | Algorithms for simulating transport and percolation        |
+----------------+------------------------------------------------------------+
| ``materials``  | A collection of predefined projects consisting of a network|
|                | with suitable topology and a geometry with necessary models|
+----------------+------------------------------------------------------------+
| ``topotools``  | Tools for querying and manipulating network topology       |
+----------------+------------------------------------------------------------+
| ``io``         | Import from and export to various common data formats      |
+----------------+------------------------------------------------------------+
| ``utils``      | Helper utilites & classes, including ``Workspace`` and     |
|                |  ``Project``                                               |
+----------------+------------------------------------------------------------+
| ``models``     | Library of pore-scale models for calculating geometric,    |
|                | thermodynamic, and physical properties                     |
+----------------+------------------------------------------------------------+

"""

__version__ = '2.1.0'

import numpy as np
np.seterr(divide='ignore', invalid='ignore')

from . import utils
from .utils import Workspace, Project
from . import core
from . import models
from . import network
from . import topotools
from . import geometry
from . import phases
from . import physics
from . import algorithms
from . import io
from . import materials
