# -*- coding: utf-8 -*-
"""
===============================================================================
Empty: Create an empty network suitable for adding to manually
===============================================================================

"""
import numpy as np
import scipy as sp
from openpnm.network import GenericNetwork


class Empty(GenericNetwork):
    r"""

    """
    def __init__(self, name=None, project=None):
        self['pore.coords'] = sp.ndarray(shape=(0, 3), dtype=float)
        self['throat.conns'] = sp.ndarray(shape=(0, 2), dtype=int)
        super().__init__(name=name, project=project)

    def _set_Np(self, Np):
        self.update({'pore.all': sp.ones(shape=(Np, ), dtype=bool)})

    def _get_Np(self):
        return sp.size(self['pore.all'])

    Np = property(fset=_set_Np, fget=_get_Np)

    def _set_Nt(self, Nt):
        self.update({'throat.all': sp.ones(shape=(Nt, ), dtype=bool)})

    def _get_Nt(self):
        return sp.size(self['throat.all'])

    Nt = property(fset=_set_Nt, fget=_get_Nt)
