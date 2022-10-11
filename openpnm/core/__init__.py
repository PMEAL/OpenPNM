r"""
Main classes of OpenPNM
=======================

This module contains the main classes from which all other major objects
(Network, Geometry, Physics, Phase, and Algorithm) derive.

The Base class
--------------

The ``Base`` class is a ``dict`` that has added methods for indexing the pores
and throats, applying labels, and managing the stored data. All OpenPNM
object inherit from ``Base`` so possess these methods.

----

``Base`` objects, Networks, Phase, Algorithms, are assigned to all locations
in the domain.  The ``Subdomain`` class is a direct descendent of ``Base``
which has the added ability to be assigned to a subset of the domain.  Objects
that inherit from ``Subdomain`` are Geomery and Physics.

Boss objects refer to the Full Domain object it is associated with.  For
Geomery objects this is the Network, and for Physics objects this is the
Phase that was specified during instantiation.

The associations between an object and it's boss are tracked using labels in
the boss.  So a Geometry object named ``geom1`` will put labels 'pore.geom1'
and 'throat.geom1' into the Network dictionary, with ``True`` values indicating
where ``geom1`` applies.

The ModelsMixin class
---------------------

`Mixins <https://en.wikipedia.org/wiki/Mixin>`_ are a useful feature of Python
that allow a few methods to be added to a class that needs them. In OpenPNM,
the ability to store and run 'pore-scale' models is not needed by some objects
(Network, Algorithms), but is essential to Geometry, Physics, and Phase
objects.

In addition to these methods, the ``ModelsMixin`` also adds a ``models``
attribute to each object.  This is a dictionary that stores the pore-scale
models and their associated parameters.  When ``regenerate_models`` is called
the function and all the given parameters are retrieved from this dictionary
and run.

"""

from ._models import *
from ._mixins import *
from ._base2 import *
