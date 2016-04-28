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
import sys as _sys
import scipy as _sp

# Check Python version
if _sys.version_info < (3, 3):
    raise Exception('OpenPNM requires Python 3.3 or greater to run')

if _sp.__version__ < '0.14.0':
    raise Exception('OpenPNM requires SciPy version 0.14.0 or greater')

__version__ = '1.4.6'

from . import Base
from . import Network
from . import Geometry
from . import Phases
from . import Physics
from . import Algorithms
from . import Postprocessing
from . import Utilities
from .Base import Controller as ctrl

_workspace = ctrl()
del ctrl
save = _workspace.save
load = _workspace.load
export_data = _workspace.export_data
import_data = _workspace.import_data
clear = _workspace.clear
purge_object = _workspace.purge_object


def view():
    print(_workspace.__str__())
