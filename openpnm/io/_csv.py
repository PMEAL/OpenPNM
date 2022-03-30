import re
import numpy as np
from openpnm.io._pandas import Pandas
from openpnm.io import GenericIO, Dict
from openpnm.utils import Workspace

ws = Workspace()


class CSV(GenericIO):
    r"""
    Reads and writes CSV (comma-separated-value files) containing pore and
    throat data

    Notes
    -----
    There are a few rules governing how the data is be stored:

    1. The first row of the file (column headers) must contain the
    property names. The subsequent rows contain the data.

    2. The property names should be in the usual OpenPNM format, such as
    of ``pore.volume`` or ``throat.surface_area``.

    3. Each column represents a specific property.  For Np x 1 or Nt x 1
    data such as *pore.volume* this is straightforward.  For Np x *m* or
    Nt x *m* data, each of the *m* columns should have their own column in
    in the CSV file, with a numpy-style index indicating which axis it
    corresponds to.  For instance, the *pore.coords* values should be stored
    as three separate columns with the headings: *pore.coords[0]*,
    *pore.coords[1]*, and *pore.coords[2]*.  OpenPNM will convert that back
    into an Np x *m* array upon loading.

    4. The file can contain both or either pore and throat data.

    5. Labels can be imported by placing the characters TRUE and FALSE
    in a column corresponding to the label name (i.e. *pore.front*).  TRUE
    indicates where the label applies and FALSE otherwise.

    """

    @classmethod
    def export_data(cls, network=None, phases=[], filename='', delim=' | '):
        r"""
        Save all the pore and throat property data on the Network (and
        optionally on any Phases objects) to CSV files.

        Parameters
        ----------
        network : OpenPNM Network
            The Network containing the data to be stored

        phases : list of OpenPNM Phases (optional)
            The Phases whose data should be stored.

        filename : str or path object
            The name of the file to store the data

        Notes
        -----
        The data from all Geometry objects is added to the file automatically.

        """
        project, network, phases = cls._parse_args(network=network,
                                                   phases=phases)
        df = Pandas.export_data(network=network, phases=phases,
                                 join=True, delim=delim)
        # Write to file
        if filename == '':
            filename = project.name
        fname = cls._parse_filename(filename=filename, ext='csv')
        df.to_csv(fname, index=False)

    @classmethod
    def import_data(cls, filename, project=None, delim=' | '):
        r"""
        Opens a 'csv' file, reads in the data, and adds it to the **Network**

        Parameters
        ----------
        filename : str (optional)
            The name of the file containing the data to import.  The formatting
            of this file is outlined below.

        project : Project
            A GenericNetwork is created and added to the specified Project.
            If no Project object is supplied then one will be created and
            returned.

        Returns
        -------
        project : list
            An OpenPNM project containing the data assigned to Generic
            versions of the objects from which it was exported.

        """
        from pandas import read_table

        if project is None:
            project = ws.new_project()

        fname = cls._parse_filename(filename, ext='csv')
        a = read_table(filepath_or_buffer=fname,
                       sep=',',
                       skipinitialspace=True,
                       index_col=False,
                       true_values=['T', 't', 'True', 'true', 'TRUE'],
                       false_values=['F', 'f', 'False', 'false', 'FALSE'])

        dct = {}
        # First parse through all the items and re-merge columns
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

        project = Dict.from_dict(dct, project=project, delim=delim)

        return project


def from_csv(filename, project=None, delim=' | '):
    project = CSV.import_data(filename=filename, project=project, delim=delim)
    return project


from_csv.__doc__ = CSV.import_data.__doc__


def to_csv(network, phases=[], filename='', delim=' | '):
    CSV.export_data(network=network, phases=phases,
                    filename=filename, delim=delim)


to_csv.__doc__ = CSV.export_data.__doc__
