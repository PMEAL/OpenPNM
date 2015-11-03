# -*- coding: utf-8 -*-
"""
===============================================================================
Import: Import networks from a standardized file format
===============================================================================

"""
import scipy as sp
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
        rarr = sp.recfromcsv(file)
        try:
            Nt = len(rarr['throat_conns'])
        except:
            raise Exception('throat_conns was not found, cannot proceed')
        try:
            Np = sp.sum([len(rarr['pore_coords'][i]) > 0 for i in range(Nt)])
        except:
            raise Exception('pore.coords was not found, cannot proceed')
        # Add basic info to Network
        self.update({'pore.all': sp.ones((Np,), dtype=bool)})
        self.update({'throat.all': sp.ones((Nt,), dtype=bool)})
        data = [sp.fromstring(rarr['pore_coords'][i], sep=' ') for i in range(Np)]
        self.update({'pore.coords': sp.vstack(data)})
        data = [sp.fromstring(rarr['throat_conns'][i], sep=' ') for i in range(Nt)]
        self.update({'throat.conns': sp.vstack(data)})
        items = list(rarr.dtype.names)
        items.remove('pore_coords')
        items.remove('throat_conns')
        # Now parse through all the other items
        for item in items:
            element = item.split('_')[0]
            prop = item.split('_')[1]
            data = rarr[item]
            if data.dtype.char == 'S':
                print('not sure how to deal with this yet')
            else:
                self.update({element+'.'+prop: data[0:self._count(element)]})
