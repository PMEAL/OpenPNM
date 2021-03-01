import numpy as np
from flatdict import FlatDict
from collections import namedtuple
from openpnm.io import Dict, GenericIO
from openpnm.utils import sanitize_dict, logging
logger = logging.getLogger(__name__)


class Pandas(GenericIO):
    r"""
    Combines all data arrays into a Pandas DataFrame object

    The structure of a DataFrame is a very close match to OpenPNMs data
    storage.  Each key becomes a column header in the Dataframe, and each
    pore or throat entry becomes a row.

    Limitations of the DataFrame are the inability to have multidimensional
    data in a single column.  The methods on a DataFrame are also oriented
    towards time-series data.

    Nonetheless, Pandas offers many useful features such as performing
    statistical analysis on property.  DataFrames also offer *many* options for
    exporting to other file formats, so if a format is not yet supported
    by OpenPNM, this could be an solution.

    """
    @classmethod
    def to_dataframe(cls, *args, **kwargs):
        r"""
        This method is being deprecated.  Use ``export_data`` instead.
        """
        data = cls.export_data(*args, **kwargs)
        return data

    @classmethod
    def export_data(cls, network=None, phases=[], join=False, delim=' | '):
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

        Returns
        -------
        Pandas ``DataFrame`` object containing property and label data in each
        column.  If ``join`` was False (default) the two DataFrames are
        returned i a named tuple, or else a single DataFrame with pore and
        throat data in the same file, despite the column length being
        different.

        """
        from pandas import DataFrame

        project, network, phases = cls._parse_args(network=network,
                                                   phases=phases)

        # Initialize pore and throat data dictionary using Dict class
        pdata = Dict.to_dict(network=network, phases=phases, element='pore',
                             interleave=True, flatten=True,
                             categorize_by=['object'])
        tdata = Dict.to_dict(network=network, phases=phases, element='throat',
                             interleave=True, flatten=True,
                             categorize_by=['object'])
        pdata = FlatDict(pdata, delimiter=delim)
        tdata = FlatDict(tdata, delimiter=delim)

        # Scan data and convert non-1d arrays to multiple columns
        for key in list(pdata.keys()):
            if np.shape(pdata[key]) != (network[0].Np,):
                arr = pdata.pop(key)
                tmp = np.split(arr, arr.shape[1], axis=1)
                cols = range(len(tmp))
                pdata.update({key+'['+str(i)+']': tmp[i].squeeze()
                              for i in cols})
        for key in list(tdata.keys()):
            if np.shape(tdata[key]) != (network[0].Nt,):
                arr = tdata.pop(key)
                tmp = np.split(arr, arr.shape[1], axis=1)
                cols = range(len(tmp))
                tdata.update({key+'['+str(i)+']': tmp[i].squeeze()
                              for i in cols})

        # Convert sanitized dictionaries to DataFrames
        pdata = DataFrame(sanitize_dict(pdata))
        tdata = DataFrame(sanitize_dict(tdata))

        # Prepare DataFrames to be returned
        if join:
            data = tdata.join(other=pdata, how='left')
        else:
            nt = namedtuple('dataframes', ('pore', 'throat'))
            data = nt(pore=pdata, throat=tdata)

        return data
