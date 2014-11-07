r"""
##################################################################################
:mod:`OpenPNM` --  A scientific pore network simulator for porous media transport
##################################################################################

Documentation is available in the docstrings and in the on-line documentation.


Subpackages
-----------

.. list-table:: OpenPNM Submodule structure.
   :widths: 10 80
   :header-rows: 1

   * - Name
     - Description
   * - :mod:`OpenPNM.Network`
     - Storage and manipulations of network topoologies and data stored on them.
   * - :mod:`OpenPNM.Geometry`
     - Geometry for pore networks. (Random cubic, image based, Voronoi). Should also contain
       a mapper of the pore network back on the segmented image.
   * - :mod:`OpenPNM.Physics`
     - Module containing pore scale physics models and equations.
   * - :mod:`OpenPNM.Phases`
     - Module containing thremophyics property estimation models.
   * - :mod:`OpenPNM.Algorithms`
     - Module containing all algorithmic classes for networks.
   * - :mod:`OpenPNM.Utilities`
     - common utilities and classes used by most of the of the modules

Import
------
>>> import OpenPNM

"""

__version__ = '1.1-beta'


__requires__ = [
    'scipy'
]

from . import Base
from . import Utilities
from . import Network
from . import Geometry
from . import Phases
from . import Physics
from . import Algorithms
from . import Postprocessing


