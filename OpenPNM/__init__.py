
r"""
##################################################################################
:mod:`OpenPNM` --  A scientific pore network simulator for porous transport media
##################################################################################
.. module:: OpenPNM
    :platform: Linux, Windows

Documentation is available in the docstrings and in ths sphinx documentation.

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
   * - :mod:`OpenPNM.Visualization`
     - Module for performing vtk-based post-processing routines.




Utility tools
-------------
::

 TODO                --- Todo


Import
------
>>> import OpenPNM as PNM

Inheritance Diagram
--------------------

.. inheritance-diagram:: OpenPNM.Network.GenericNetwork

Package Documentation
---------------------

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

.. automodule:: Visualization
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: GUI
   :members:
   :undoc-members:
   :show-inheritance:
"""

__version__ = '1.0.0'

__requires__ = [
    'scipy'
]

from .__Base__ import Base
from .__Core__ import Core
from . import Utilities
from . import Network
from . import Geometry
from . import Phases
from . import Physics
from . import Algorithms
from . import Postprocessing
from .Utilities.__Load__ import Load





