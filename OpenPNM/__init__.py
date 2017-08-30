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

Example
-------
The following few lines setup all the necessary objects for performing
simulations:

>>> import OpenPNM as op
>>> net = op.Network.Cubic(shape=[5, 5, 5], spacing=0.001)
>>> geo = op.Geometry.Stick_and_Ball(network=net, pores=net.Ps, throats=net.Ts)
>>> air = op.Phases.Air(network=net)
>>> water = op.Phases.Water(network=net)
>>> phys_air = op.Physics.Standard(network=net, geometry=geo, phase=air)
>>> phys_wat = op.Physics.Standard(network=net, geometry=geo, phase=water)

"""
import sys as _sys
import scipy as _sp

# Check Python version
if _sys.version_info < (3, 4):
    raise Exception('OpenPNM requires Python 3.4 or greater to run')

__version__ = '1.6.2'

from . import Base
from . import Network
from . import Geometry
from . import Phases
from . import Physics
from . import Algorithms
from . import Postprocessing
from . import Utilities
from .Base import Workspace as mgr

_workspace = mgr()
del mgr
save_workspace = _workspace.save_workspace
load_workspace = _workspace.load_workspace
export_data = _workspace.export_data
import_data = _workspace.import_data
clear = _workspace.clear
purge_object = _workspace.purge_object


def view():
    print(_workspace.__str__())
