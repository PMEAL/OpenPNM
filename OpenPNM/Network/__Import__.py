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
            N = self._count(element)
            prop = item.split('_')[1]
            data = rarr[item]
            if data.dtype.char == 'S':
                if data[0].decode().upper()[0] in ['T', 'F']: # If data is True and False
                    data = sp.chararray.decode(data)
                    data = sp.chararray.upper(data)
                    ind = sp.where(data == 'T')[0]
                    data = sp.zeros((N,), dtype=bool)
                    data[ind] = True
                else:  # If data is an array of lists
                    data = [list(sp.fromstring(rarr[item][i], sep=' ')) for i in range(N)]
                    data = sp.array(data)
            self.update({element+'.'+prop: data[0:N]})
