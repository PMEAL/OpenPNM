r"""

**openpnm.materials**

----

This module provides a library of preconfigured Network-Geometry combinations.

In most case the topology and geometry cannot be considered in isolation.
This module provides recipes that create both the Network and Geometry objects
simultaneously to ensure sensible correspondance between things like lattice
spacing and pore sizes.  Some of the classes in this module have a signficant
amount of custom code (e.g. ``VoronoiFibers``), while others are simple
recipes that combine existing models in OpenPNM (e.g. ``CubicSandstone``).


The table below gives a list of available Material generators:

+---------------------+-------------------------------------------------------+
| Material Name       | Description                                           |
+=====================+=======================================================+
| VoronoiFibers       | Resembles a fibrous paper or mat with straight        |
|                     | intersecting fibers.                                  |
+---------------------+-------------------------------------------------------+
| BundleOfTubes       | A 2D network with cylindrical pores through the       |
|                     | thickness                                             |
+---------------------+-------------------------------------------------------+
| CubicSandstone      | A cubic 3D network representing common sandstones     |
+---------------------+-------------------------------------------------------+

"""

from .VoronoiFibers import VoronoiFibers
from .BundleOfTubes import BundleOfTubes
from .CubicSandstone import CubicSandstone
