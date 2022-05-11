import numpy as np
from pandas import DataFrame
from openpnm.io import Dict
from openpnm.utils import sanitize_dict


def project_to_pandas(project, join=False, delim='.'):
    r"""
    Convert the Network and Phase data to a Pandas DataFrame(s)

    Parameters
    ----------
    project : list
        An OpenPNM ``project`` object
    join : bool
        If ``False`` (default), two DataFrames are returned with *pore*
        data in one, and *throat* data in the other.  If ``True`` the pore
        and throat data are combined into a single DataFrame.  This can be
        problematic as it will put NaNs into all the *pore* columns which
        are shorter than the *throat* columns.

    Returns
    -------
    Pandas ``DataFrame`` object(s) containing property and label data in
    each column. If ``join`` was ``False`` (default) the two DataFrames are
    returned in a dict, otherwise a single DataFrame with pore and
    throat data in the same file, despite the column length being
    different.

    """
    network = project.network
    phases = project.phases

    # Initialize pore and throat data dictionary using Dict class
    pdata = Dict.to_dict(network=network, phases=phases, element='pore',
                         flatten=True, categorize_by=[], delim=delim)
    tdata = Dict.to_dict(network=network, phases=phases, element='throat',
                         flatten=True, categorize_by=[], delim=delim)

    # Scan data and convert non-1d arrays to multiple columns
    for key in list(pdata.keys()):
        if np.shape(pdata[key]) != (network.Np,):
            arr = pdata.pop(key)
            tmp = np.split(arr, arr.shape[1], axis=1)
            cols = range(len(tmp))
            pdata.update({key+'['+str(i)+']': tmp[i].squeeze()
                          for i in cols})
    for key in list(tdata.keys()):
        if np.shape(tdata[key]) != (network.Nt,):
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
        data = {'pore': pdata, 'throat': tdata}
    return data
