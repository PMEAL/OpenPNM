r"""
Main classes of OpenPNM
=======================

This module contains the main classes from which all other major objects derive.

The Base class
--------------

The ``Base`` class is a ``dict`` that has added methods for indexing the pores
and throats, applying labels, and managing the stored data. All OpenPNM
object inherit from ``Base`` so possess these methods.

"""

from ._models import *
from ._mixins import *
from ._base2 import *
