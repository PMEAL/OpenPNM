# -*- coding: utf-8 -*-
"""
===============================================================================
Import: Import networks from a standardized file format
===============================================================================

"""
import scipy as sp
from OpenPNM.Network import GenericNetwork
from OpenPNM.Utilities import topology
from OpenPNM.Utilities import IO as io
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)
topo = topology()


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
