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
    or Pandas

    2. **MAT** : Matlab files containing named arrays

    3. **Pajek** : A native output format of NetworkX

    4. **VTK** : Visualization Tool-kit files in the legacy format.

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def from_csv(self, filename, overwrite=True):
        io.CSV.load(network=self, filename=filename, overwrite=overwrite)

    from_csv.__doc__ = io.CSV.load.__doc__

    def from_vtk(self, filename):
        io.VTK.load(filename=filename)

    def from_mat(self, filename, overwrite=True):
        io.MAT.load(network=self, filename=filename, overwrite=overwrite)

    from_mat.__doc__ = io.MAT.load.__doc__

    def from_pajek(self, filename, overwrite=True):
        raise NotImplemented
        # io.PAJEK.load(network=self, filename=filename, overwrite=overwrite)
