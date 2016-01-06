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
    This class is used to import data resulting from some external network
    extraction tools {refs here?}.  The aim of this class is to define a
    standard way to represent network data and transfer into OpenPNM.   The
    main principle is to keep is as simple and general as possible, so the
    CSV data format was used.

    There are a few rules governing how the data should be stored:

    1.  The first row of the file (column headers) should contain the property
    name

    2.  The property names should be in the format of *pore_volume* or
    *throat_length*.  In OpenPNM this will become *pore.volume* or
    *throat.length* (i.e. the underscore is replaced by a dot).

    3.  Each column represents a specific property.  For Np x 1 or Nt x 1 data
    such as *pore_volume* this is straightforward.  For Np x m or Nt x m data,
    it must be entered in as a set of values NOT separated by commas.  For
    instance, the *pore_coords* values should be X Y Z with spaces, not commas
    between them.

    4.  Labels can also be imported by placing the characters T and F in a
    column corresponding to the label name (i.e. *pore_front*).  T indicates
    where the label applies and F otherwise.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def from_csv(self, file):
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
        data = [sp.fromstring(rarr['pore_coords'][i], sep=' ')
                for i in range(Np)]
        self.update({'pore.coords': sp.vstack(data)})
        data = [sp.fromstring(rarr['throat_conns'][i], sep=' ')
                for i in range(Nt)]
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
                if data[0].decode().upper()[0] in ['T', 'F']:  # If boolean
                    data = sp.chararray.decode(data)
                    data = sp.chararray.upper(data)
                    ind = sp.where(data == 'T')[0]
                    data = sp.zeros((N,), dtype=bool)
                    data[ind] = True
                else:  # If data is an array of lists
                    data = [list(sp.fromstring(rarr[item][i], sep=' '))
                            for i in range(N)]
                    data = sp.array(data)
            self.update({element+'.'+prop: data[0:N]})
