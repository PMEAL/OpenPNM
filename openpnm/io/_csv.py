import re

import numpy as np
from pandas import read_table

from openpnm.io import _parse_filename
from openpnm.io._pandas import network_to_pandas, project_to_pandas
from openpnm.utils import Project


def project_to_csv(project, filename=''):
    r"""
    Save all the pore and throat data on the Network and Phase objects to a CSV
    file

    Parameters
    ----------
    project : list
        An openpnm ``project`` object
    filename : str or path object
        The name of the file to store the data

    """
    df = project_to_pandas(project=project, join=True, delim='.')
    if filename == '':
        filename = project.name
    fname = _parse_filename(filename=filename, ext='csv')
    df.to_csv(fname, index=False)


def network_to_csv(network, filename=''):
    """Exports a network to a CSV file."""
    proj = Project()
    proj.append(network)
    df = network_to_pandas(network=network, join=True, delim='.')
    if filename == '':
        filename = network.name
    fname = _parse_filename(filename=filename, ext='csv')
    df.to_csv(fname, index=False)


def network_from_csv(filename):
    """Loads a network from a CSV file."""
    from openpnm.network import Network
    fname = _parse_filename(filename=filename, ext='csv')

    a = read_table(filepath_or_buffer=fname,
                   sep=',',
                   skipinitialspace=True,
                   index_col=False,
                   true_values=['T', 't', 'True', 'true', 'TRUE'],
                   false_values=['F', 'f', 'False', 'false', 'FALSE'])

    # First parse through all the items and re-merge columns
    dct = {}
    keys = sorted(list(a.keys()))
    for item in keys:
        m = re.search(r'\[.\]', item)  # The dot '.' is a wildcard
        if m:  # m is None if pattern not found, otherwise merge cols
            pname = re.split(r'\[.\]', item)[0]  # Get base propname
            # Find all other keys with same base propname
            merge_keys = [k for k in a.keys() if k.startswith(pname)]
            # Rerieve and remove arrays with same base propname
            merge_cols = [a.pop(k) for k in merge_keys]
            # Merge arrays into multi-column array and store in DataFrame
            dct[pname] = np.vstack(merge_cols).T
            # Remove key from list of keys
            for k in keys:
                if k.startswith(pname):
                    keys.pop(keys.index(k))
        else:
            dct[item] = np.array(a.pop(item))

    # Now scan through 'pore.coords' and 'throat.conns' to get Np and Nt,
    # then remove the nans
    try:
        Np = np.where(np.isnan(dct['pore.coords'][:, 0]))[0][0]
    except IndexError:
        Np = dct['pore.coords'][:, 0].shape[0]
    try:
        Nt = np.where(np.isnan(dct['throat.conns'][:, 0]))[0][0]
    except IndexError:
        Nt = dct['throat.conns'][:, 0].shape[0]
    for k, v in dct.items():
        if k.startswith('pore.'):
            dct[k] = v[:Np, ...]
        if k.startswith('throat.'):
            dct[k] = v[:Nt, ...]

    network = Network()
    network.update(dct)
    return network
