import logging
from openpnm.io import project_to_dict, _parse_filename
import pandas as pd


logger = logging.getLogger(__name__)


def project_to_hdf5(project, filename=''):
    r"""
    Creates an HDF5 file containing data from the specified objects

    Parameters
    ----------
    network : Network
        The network containing the desired data

    phases : list[Phase]s (optional, default is none)
        A list of phase objects whose data are to be included

    Returns
    -------
    f : hdf5 file handle
        A handle to an hdf5 file.  This must be closed when done (i.e.
        ``f.close()``.
    """
    from h5py import File as hdfFile
    if filename == '':
        filename = project.name
    filename = _parse_filename(filename, ext='hdf')

    dct = project_to_dict(project=project)
    d = pd.json_normalize(dct, sep='.').to_dict(orient='records')[0]
    for k in list(d.keys()):
        d[k.replace('.', '/')] = d.pop(k)
    f = hdfFile(filename, "w")
    for item in list(d.keys()):
        tempname = '_'.join(item.split('.'))
        arr = d[item]
        if d[item].dtype == 'O':
            logger.warning(item + ' has dtype object, will not write to file')
            del d[item]
        elif 'U' in str(arr[0].dtype):
            pass
        else:
            f.create_dataset(name='/'+tempname, shape=arr.shape,
                             dtype=arr.dtype, data=arr)
    return f


def print_hdf5(f, flat=False):
    r"""
    Given an hdf5 file handle, prints to console in a human readable manner

    Parameters
    ----------
    f : file handle
        The hdf5 file to print
    flat : bool
        Flag to indicate if print should be nested or flat.  The default is
        ``flat==False`` resulting in a nested view.
    """
    if flat is False:
        def print_level(f, p='', indent='-'):
            for item in f.keys():
                if hasattr(f[item], 'keys'):
                    p = print_level(f[item], p=p, indent=indent + indent[0])
                elif indent[-1] != ' ':
                    indent = indent + ''
                p = indent + item + '\n' + p
            return p
        p = print_level(f)
        print(p)
    else:
        f.visit(print)
