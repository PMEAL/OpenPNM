r"""
===
I/O
===

This module contains functionality for importing and exporting data
between OpenPNM and other formats.

Note that saving and loading OpenPNM data in its own format is done
through the Workspace and Project objects.  Saving data through these
classes losses lots of vital information about the simulation, such as
pore scale-models and parameters.

"""
from ._generic import GenericIO
from ._pickle import Pickle
from ._dict import Dict
from ._vtk import VTK
from ._csv import CSV
from ._networkx import NetworkX
from ._mat import MAT
from ._imorph import iMorph
from ._marock import MARock
from ._statoil import Statoil
from ._pergeos import PerGeos
from ._porespy import PoreSpy
from ._pandas import Pandas
from ._hdf5 import HDF5
from ._xdmf import XDMF
from ._jsongraph import JSONGraphFormat
from ._stl import STL
from ._comsol import COMSOL
from ._salome import Salome
from ._pnm import PNM
from ._paraview import ParaView
