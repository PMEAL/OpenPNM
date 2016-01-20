# -*- coding: utf-8 -*-
"""
===============================================================================
Import: Import networks from a standardized file format
===============================================================================

"""
from OpenPNM.Network import GenericNetwork
from OpenPNM.Utilities import IO as io
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class Import(GenericNetwork):
    r"""
    This method is used to import data resulting from some external network
    extraction tools onto an OpenPNM Network object.

    Notes
    -----
    This class currently has support for:

    1. **CSV** : Comma-separated values such as spreadsheets saved from Excel
    or `Pandas <http://pandas.pydata.org/>`_

    2. **MAT** : Matlab `'Mat-files'
    <http://www.mathworks.com/help/matlab/ref/matfile.html>`_ containing
    named arrays

    3. **YAML** : An output format of `NetworkX <https://networkx.github.io/>`_

    4. **VTK** : `Visualization Toolkit <http://www.vtk.org/>`_ files in the
    *vtp* or unstructure grid format.

    """

    def from_csv(self, filename, mode='overwrite'):
        io.CSV.load(network=self, filename=filename, mode=mode)

    def from_vtk(self, filename, mode='overwrite'):
        io.VTK.load(network=self, filename=filename, mode=mode)

    def from_mat(self, filename, mode='overwrite'):
        io.MAT.load(network=self, filename=filename, mode=mode)

    def from_yaml(self, filename, mode='overwrite'):
        io.YAML.load(network=self, filename=filename, mode=mode)
