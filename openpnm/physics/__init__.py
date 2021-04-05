r"""

**openpnm.physics**

----

This module contains the GenericPhysics class, along with a few preconfigured
classes that include common pore-scale physical models

----

**The GenericPhysics Class**

----

**Library of Preconfigured Physics Classes**

+-------------+---------------------------------------------------------------+
| Class       | Comments                                                      |
+=============+===============================================================+
| Standard    |                                                               |
+-------------+---------------------------------------------------------------+

"""

from .GenericPhysics import GenericPhysics
from .Standard import Standard
from .Standard2D import Standard2D
from .Classic import Classic
from .Basic import Basic
