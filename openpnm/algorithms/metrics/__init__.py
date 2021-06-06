r"""

**openpnm.algorithms.metrics**

----

This module provides a library of preconfigured algorithms.


The table below gives a list of available Metrics:

+---------------------+-------------------------------------------------------+
| Material Name       | Description                                           |
+=====================+=======================================================+
| MercuryIntrusion    | Invades mercury into a network from all sides         |
+---------------------+-------------------------------------------------------+
| FormationFactor     | Calculates the effective conductivity of a network,   |
|                     | and normalizes it by the solution conductivity        |
+---------------------+-------------------------------------------------------+
| MercuryIntrusion    | Invades mercury into a network from all sides         |
+---------------------+-------------------------------------------------------+

"""

from .GenericMetric import GenericMetric
from .MercuryIntrusion import MercuryIntrusion
from .FormationFactor import FormationFactor
from .RelativePermeability import RelativePermeability
from .PNFlow import PNFlow
