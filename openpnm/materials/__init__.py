r"""

**openpnm.materials**

----

This module provides a library of preconfigured Network-Geometry combinations.

In most case the topology and geometry cannot be considered in isolation.
This module provides recipes that create both the Network and Geometry objects
simultaneously to ensure sensible correspondance between things like lattice
spacing and pore sizes.  Some of the classes in this module have a signficant
amount of custom code (e.g. ``VoronoiFibers``), while others are simple
recipes that combine existing models in OpenPNM (e.g. ``BereaCubic``).


The table below gives a list of available Material generators:

+---------------------+-------------------------------------------------------+
| Material Name       | Description                                           |
+=====================+=======================================================+
| CubicBerea          | A traditional Berea Sandstone on a Cubic lattice      |
+---------------------+-------------------------------------------------------+
| VoronoiFibers       | Resembles a fibrous paper or mat with straight        |
|                     | intersecting fibers.                                  |
+---------------------+-------------------------------------------------------+

Examples
--------

>>> import openpnm as op
>>> proj = op.materials.BereaCubic(shape=[5, 5, 5])

The Materials classes create a Project instance, so the Network and Geometry
objects must be retrieved explicitly:

>>> net = proj.network
>>> print(net.props())
――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
1       : pore.coords
2       : throat.conns
――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
>>> print(proj.geometries().keys())
dict_keys(['geo_01'])
>>> geo = proj.geometries()['geo_01']
>>> print(geo)
――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
openpnm.geometry.GenericGeometry:       geo_01
――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
#     Properties                                    Valid Values
――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
1     pore.diameter                                   125 / 125
2     pore.seed                                       125 / 125
3     pore.volume                                     125 / 125
4     throat.diameter                                1036 / 1036
5     throat.length                                  1036 / 1036
6     throat.volume                                  1036 / 1036
――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
#     Labels                                        Assigned Locations
――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
1     pore.all                                      125
2     throat.all                                    1036
――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――


"""
from .VoronoiFibers import VoronoiFibers
from .BereaCubic import BereaCubic
