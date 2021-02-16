r"""

**openpnm.io**

----

This module contains functionality for importing and exporting data between
OpenPNM and other formats.

Note that saving and loading OpenPNM data in its own format is done through
the Workspace and Project objects.  Saving data through these classes losses
lots of vital information about the simulation, such as pore scale-models and
parameters.

----

**Available Import and Export Options**

+----------+------------------------------------------------------------------+
| Format   | Description                                                      |
+==========+==================================================================+
| Dict     | Generates hierarchical ``dicts`` with a high degree of control   |
|          | over the structure                                               |
+----------+------------------------------------------------------------------+
| Pandas   | Combines all data arrays into a Pandas DataFrame object          |
+----------+------------------------------------------------------------------+
| Pickle   | Read and write OpenPNM Workspaces, Projects, objects as Pickles  |
+----------+------------------------------------------------------------------+
| CSV      | Reads and writes CSV (comma-separated-value files) containing    |
|          | pore and throat data                                             |
+----------+------------------------------------------------------------------+
| HDF5     | The HDF5 (Hierarchical Data Format) file is good for high-       |
|          | peformance, long term data storage                               |
+----------+------------------------------------------------------------------+
| XDMF     | The eXtensible Data Model Format combines XML descriptors with   |
|          | HDF5 data storage                                                |
+----------+------------------------------------------------------------------+
| VTK      | The Visualization Toolkit (VTK) format defined by Kitware and    |
|          | used by Paraview                                                 |
+----------+------------------------------------------------------------------+
| NetworkX | NetworkX is a common tool for dealing with network structures    |
+----------+------------------------------------------------------------------+
| MAT      | MAT files are a format used by Matlab                            |
+----------+------------------------------------------------------------------+
| JGF      | The JSON Graph Format is a schema specification of how to write  |
|          | a graph object into JSON format                                  |
+----------+------------------------------------------------------------------+
| STL      | The STL Format is a Standard Triangle (or Tessellation) Language |
|          | supported by many CAD packages and used for 3D printing          |
+----------+------------------------------------------------------------------+
| Salome   | A .py instruction file to be read by Salome to build a geomertry |
+----------+------------------------------------------------------------------+
| PerGeos  | The PerGeos format is used by the commercial software Avizo      |
+----------+------------------------------------------------------------------+
| PoreSpy  | PoreSpy contains the snow extraction algorithm                   |
+----------+------------------------------------------------------------------+
| iMorph   | iMorph is a graphical interface program that provides some image |
|          | analysis tools for porous media                                  |
+----------+------------------------------------------------------------------+
| MARock   | 3DMA-Rock is a network extraction algorithm                      |
+----------+------------------------------------------------------------------+
| Statoil  | The StatOil format is used by the Maximal Ball network extraction|
|          | code of the Imperial College London group                        |
+----------+------------------------------------------------------------------+

Of the above, the last three only support *reading* from files and there are
no plans to implement writing such files.

"""
from .GenericIO import GenericIO
from .Pickle import Pickle
from .Dict import Dict
from .VTK import VTK
from .CSV import CSV
from .NetworkX import NetworkX
from .MAT import MAT
from .iMorph import iMorph
from .MARock import MARock
from .Statoil import Statoil
from .PerGeos import PerGeos
from .PoreSpy import PoreSpy
from .Pandas import Pandas
from .HDF5 import HDF5
from .XDMF import XDMF
from .JSONGraphFormat import JSONGraphFormat
from .STL import STL
from .COMSOL import COMSOL
from .Salome import Salome
from .PNM import PNM
from .ParaView import ParaView
