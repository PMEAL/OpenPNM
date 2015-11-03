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
        if type(file) is str:
            pass
        rec_arr = sp.recfromcsv(file)
        Np = sp.sum(~sp.isnan(rec_arr['coords1']))
        Nt = sp.sum(~sp.isnan(rec_arr['conns1']))
        self.update({'pore.all': sp.ones((Np,), dtype=bool)})
        self.update({'throat.all': sp.ones((Nt,), dtype=bool)})
        for item in rec_arr.dtype.names:
            data = rec_arr[item]
            if sp.sum(~sp.isnan(data)) == Np:
                self.update({'pore.'+item: data})
            if sp.sum(~sp.isnan(data)) == Nt:
                self.update({'throat.'+item: data})
