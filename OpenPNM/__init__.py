
r"""
##################################################################################
:mod:`OpenPNM` --  A scientific pore network simulator for porous media transport
##################################################################################
.. module:: OpenPNM
    :platform: Linux, Windows

Documentation is available in the docstrings and in the sphinx documentation.

Contents
--------
The OpenPNM package imports all the functions from the top level modules.


Subpackages
-----------


.. list-table:: OpenPNM submodule structure.
   :widths: 10 80
   :header-rows: 1

   * - Name
     - Description
   * - :mod:`OpenPNM.Utilities`
     - common utilities and classes used by most of the of the modules
   * - :mod:`OpenPNM.Network`
     - Storage and manipulations of network topoologies and data stored on them.
   * - :mod:`OpenPNM.Geometry`
     - Geometry for pore networks. (Random cubic, image based, Voronoi). Should also contain
       a mapper of the pore network back on the segmented image.
   * - :mod:`OpenPNM.Algorithms`
     - Module containing all algorithmic classes for networks.
   * - :mod:`OpenPNM.Physics`
     - Module containing pore scale physics models and equations.


Import
------
>>> import OpenPNM as PNM

Inheritance Diagram
--------------------

.. inheritance-diagram:: OpenPNM.Network.GenericNetwork

Package Documentation
---------------------

.. automodule:: Base
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Utilities
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Network
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Geometry
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Algorithms
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Physics
   :members:
   :undoc-members:
   :show-inheritance:
   
.. automodule:: Phases
   :members:
   :undoc-members:
   :show-inheritance:

"""

#__version__ = '1.0.1'

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






from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
