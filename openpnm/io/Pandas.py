import scipy as _sp
from collections import namedtuple
import pandas as _pd
from openpnm.core import logging
from openpnm.io import Dict
from openpnm.utils import FlatDict
logger = logging.getLogger(__name__)


class Pandas():

    @classmethod
    def to_dataframe(cls, network=None, phases=[], join=False):
        r"""
        Convert the Network (and optionally Phase) data to Pandas DataFrames.

        Parameters
        ----------
        network: OpenPNM Network Object
            The network containing the data to be stored

        phases : list of OpenPNM Phase Objects
            The data on each supplied phase will be added to DataFrame

        join : boolean
            If ``False`` (default), two DataFrames are returned with *pore*
            data in one, and *throat* data in the other.  If ``True`` the pore
            and throat data are combined into a single DataFrame.  This can be
            problematic as it will put NaNs into all the *pore* columns which
            are shorter than the *throat* columns.

        """
        simulation, network, phases = cls._parse_args(network=network,
                                                      phases=phases)

        # Initialize pore and throat data dictionary using Dict class
        pdata = Dict.to_dict(network=network, phases=phases, element='pore',
                             interleave=True, flatten=True,
                             categorize_by=['object'])
        tdata = Dict.to_dict(network=network, phases=phases, element='throat',
                             interleave=True, flatten=True,
                             categorize_by=['object'])
        pdata = FlatDict(pdata, delimiter=' | ')
        tdata = FlatDict(tdata, delimiter=' | ')

        # Scan data and convert non-1d arrays to multiple columns
        for key in list(pdata.keys()):
            if _sp.shape(pdata[key]) != (network.Np,):
                arr = pdata.pop(key)
                tmp = _sp.split(arr, arr.shape[1], axis=1)
                cols = range(len(tmp))
                pdata.update({key+'['+str(i)+']': tmp[i].squeeze()
                              for i in cols})
        for key in list(tdata.keys()):
            if _sp.shape(tdata[key]) != (network.Nt,):
                arr = tdata.pop(key)
                tmp = _sp.split(arr, arr.shape[1], axis=1)
                cols = range(len(tmp))
                tdata.update({key+'['+str(i)+']': tmp[i].squeeze()
                              for i in cols})

        # Convert sanitized dictionaries to DataFrames
        pdata = _pd.DataFrame(pdata)
        tdata = _pd.DataFrame(tdata)

        # Prepare DataFrames to be returned
        if join:
            data = tdata.join(other=pdata, how='left')
        else:
            data = namedtuple('dataframes', ('pore', 'throat'))
            data(pore=pdata, throat=tdata)

        return data
