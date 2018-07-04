r"""

**openpnm.core**

----

This module contains the main classes from which all other major objects
(Network, Geometry, Physics, Phase, and Algorithm) derive.

----

**The Base Class**

The ``Base`` class is a ``dict`` that has added methods for indexing the pores
and throats, applying labels, and managing the stored data.  All OpenPNM
object inherit from ``Base`` so possess these methods:

+----------------------+-------------------------------------------------+
| Method               | Description                                     |
+======================+=================================================+
| num_pores,           | Returns the number of pores and throats with the|
| num_throats          | specified labels                                |
+----------------------+-------------------------------------------------+
| props, labels        | Returns a list of all property (numeric) or     |
|                      | label (boolean) arrays on the object            |
+----------------------+-------------------------------------------------+
| pores, throats       | Returns pore or throat indicies where given     |
|                      | labels have been applied                        |
+----------------------+-------------------------------------------------+
| toindices, tomaks    | Convert a boolean mask to a list of indices and |
|                      | vice-versa                                      |
+----------------------+-------------------------------------------------+
| filter_by_label      | Reduces a list of indices based on labels       |
+----------------------+-------------------------------------------------+
| interpolate_data     | Determines a pore (or throat) property as the   |
|                      | average if it's neighbors' values               |
+----------------------+-------------------------------------------------+
| map_pores,           | Translates pore and throat ids into indices     |
| map_throats          |                                                 |
+----------------------+-------------------------------------------------+
| show_hist            | Show a quick plot of key property distributi... |
+----------------------+-------------------------------------------------+
| check_data_health    | Check the health of pore and throat data        |
+----------------------+-------------------------------------------------+

----

**The Subdomain Class**

``Base`` objects, Networks, Phases, Algorithms, are assigned to all locations
in the domain.  The ``Subdomain`` class is a direct descendent of ``Base``
which has the added ability to be assigned to a subset of the domain.  Objects
that inherit from ``Subdomain`` are Geomery and Physics.

Only two methods are added to the ``Subdomain`` class to make this work:

+------------------+-----------------------------------------------------+
| Method           | Description                                         |
+==================+=====================================================+
| drop_locations   | Removes association between an object and it's boss |
+------------------+-----------------------------------------------------+
| add_locations    | Adds associations between an object and it's boss   |
+------------------+-----------------------------------------------------+

Boss objects refer to the Full Domain object it is associated with.  For
Geomery objects this is the Network, and for Physics objects this is the
Phase that was specified during instantiation.

The associations between an object and it's boss are tracked using labels in
the boss.  So a Geometry object named ``geom1`` will put labels 'pore.geom1'
and 'throat.geom1' into the Network dictionary, with ``True`` values indicating
where ``geom1`` applies.

----

**The ModelsMixin Class**

`Mixins <https://en.wikipedia.org/wiki/Mixin>`_ are a useful feature of Python
that allow a few methods to be added to a class that needs them.  In OpenPNM,
the ability to store and run 'pore-scale' models is not needed by some objects
(Network, Algorithms), but is essential to Geometry, Physics, and Phase
objects.

The ``ModelsMixin`` adds the following few methods:

+----------------------+-------------------------------------------------+
| Method               | Description                                     |
+======================+=================================================+
| add_model            | Adds a new model to the models dictionary       |
+----------------------+-------------------------------------------------+
| remove_model         | Removes model and data from object              |
+----------------------+-------------------------------------------------+
| regenerate_models    | Re-runs the specified model or models           |
+----------------------+-------------------------------------------------+

In addition to these methods, the ``ModelsMixin`` also adds a ``models``
attribute to each object.  This is a dictionary that stores the pore-scale
models and their associated parameters.  When ``regenerate_models`` is called
the function and all the given parameters are retrieved from this dictionary
and run.


"""

from .ModelsMixin import ModelsMixin, ModelsDict
from .Base import Base
from .Subdomain import Subdomain
