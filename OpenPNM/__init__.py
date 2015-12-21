r"""
###############################################################################
:mod:`OpenPNM` --  A pore network modeling framework for simulating tranport in
porous media
###############################################################################

Documentation is available in the docstrings and in the on-line documentation.


Subpackages
-----------

.. list-table:: OpenPNM Submodule structure
   :widths: 10 80
   :header-rows: 1

   * - Name
     - Description
   * - :mod:`OpenPNM.Network`
     - Storage of network topoology and methods for querying topology.
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
import sys
import scipy as sp

# Check Python version
if sys.version_info < (3, 3):
    raise Exception('OpenPNM requires Python 3.3 or greater to run')

if sp.__version__ < '0.14.0':
    raise Exception('OpenPNM requires SciPy version 0.14.0 or greater')

__requires__ = ['scipy']
__version__ = '1.3.0'

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
clear = _controller.clear
purge_object = _controller.purge_object


def view():
    print(_controller.__str__())
