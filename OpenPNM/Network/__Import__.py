# -*- coding: utf-8 -*-
"""
===============================================================================
Import: Import networks from a standardized file format
===============================================================================

"""
import scipy as sp
import OpenPNM.Utilities.misc as misc
from OpenPNM.Network import GenericNetwork
from OpenPNM.Utilities import topology
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)
topo = topology()


class Import(GenericNetwork):
    r"""
    """
    def __init__(self, file, **kwargs):
        super().__init__(**kwargs)
