r"""
##################################################################################
:mod:`OpenPNM` --  A pore network modeling framework for simulating tranport in porous media
##################################################################################

Documentation is available in the docstrings and in the on-line documentation.


Subpackages
-----------

.. list-table:: OpenPNM Submodule structure
   :widths: 10 80
   :header-rows: 1

   * - Name
     - Description
   * - :mod:`OpenPNM.Network`
     - Storage of network topoology and methods for working with the data stored on them.
   * - :mod:`OpenPNM.Geometry`
     - Manages the geometrical properties of the pores and throats.
   * - :mod:`OpenPNM.Physics`
     - Module containing pore scale physics models and equations.
   * - :mod:`OpenPNM.Phases`
     - Module containing thremophyical property estimation models.
   * - :mod:`OpenPNM.Algorithms`
     - Module containing algorithms for performing simulations on networks.
   * - :mod:`OpenPNM.Utilities`
     - Common utilities and classes used by most of the of the modules.

Import
------
>>> import OpenPNM

"""

import scipy as sp

if sp.__version__ < '0.14.0':
    raise Exception('OpenPNM requires SciPy version 0.14.0 or greater')

__version__ = '1.2.0'

__requires__ = ['scipy']

from . import Base
from . import Utilities
from . import Network
from . import Geometry
from . import Phases
from . import Physics
from . import Algorithms
from . import Postprocessing
from .Base import Controller as ctrl

_controller = ctrl()
del ctrl
save = _controller.save
load = _controller.load
export = _controller.export
view = _controller.__str__
