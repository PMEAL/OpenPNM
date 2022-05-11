r"""
Collection of functions for import/export-ing data
==================================================

This module contains functionality for importing and exporting data
between OpenPNM and other formats.

Note that saving and loading OpenPNM data in its own format is done
through the Workspace and Project objects.  Saving data through these
classes losses lots of vital information about the simulation, such as
pore scale-models and parameters.

"""

from ._utils import *
from ._generic_io import GenericIO
from ._dict import Dict
from ._vtk import project_to_vtk
from ._pandas import project_to_pandas
from ._csv import project_to_csv
from ._hdf5 import project_to_hdf5, print_hdf5
from ._xdmf import project_to_xdmf
from ._marock import network_from_marock
from ._porespy import network_from_porespy
from ._networkx import network_to_networkx, network_from_networkx
from ._statoil import network_to_statoil, network_from_statoil
from ._pergeos import project_to_pergeos, network_from_pergeos
from ._jsongraph import network_to_jsongraph, network_from_jsongraph
# from ._stl import to_stl
# from ._comsol import to_comsol
# from ._salome import to_salome
