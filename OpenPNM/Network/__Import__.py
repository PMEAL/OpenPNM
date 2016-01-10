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
    extraction tools.  The aim of this class is to define a standard way to
    represent network data and transfer into OpenPNM.   The main principle is
    to keep it as simple and general as possible, so the CSV data format was
    used.

    """

    def __init__(self, filename=None, **kwargs):
        super().__init__(**kwargs)
        if filename:
            self.from_csv(filename=filename)

    def from_csv(self, filename, overwrite=True):
        io.CSV.load(network=self, filename=filename, overwrite=overwrite)

    from_csv.__doc__ = io.CSV.load.__doc__

    def from_vtk(self, filename):
        io.VTK.load(filename=filename)

    def from_mat(self, filename, overwrite):
        io.MAT.load(network=self, filename=filename, overwrite=overwrite)

    from_mat.__doc__ = io.MAT.load.__doc__
