import re
import scipy as sp
import pandas as pd
from openpnm.core import logging, Project
from openpnm.network import GenericNetwork
from openpnm.io import GenericIO
from openpnm.io.Pandas import Pandas
logger = logging.getLogger(__name__)


class CSV(GenericIO):
    r"""
    This class is used for reading and writing CSV files containing pore and
    throat property data.  This class uses Pandas for transferring data from
    the OpenPNM format to CSV.

    Notes
    -----
    There are a few rules governing how the data should be stored:

    1. The first row of the file (column headers) must contain the
    property names. The subsequent rows contain the data.

    2. The property names should be in the usual OpenPNM format, such as
    of *pore.volume* or *throat.surface_area*.

    3. It's possible to sub-categorize each property by prefixing the
    category to the property name, separated by forward slashes (/).  Ex:
    ``'pore.diameter'`` can be associated with a certain network using
    ``network_1/pore.diameter``.

    3. Each column represents a specific property.  For Np x 1 or Nt x 1
    data such as *pore.volume* this is straightforward.  For Np x *m* or
    Nt x *m* data, it must be entered in as a set of values NOT separated by
    commas.  For instance, the *pore.coords* values should be X Y Z with
    *spaces* between each value, not commas.

    4. The file can contain both or either pore and throat data.

    5. Labels can be imported by placing the characters TRUE and FALSE
    in a column corresponding to the label name (i.e. *pore.front*).  TRUE
    indicates where the label applies and FALSE otherwise.

    """

    @classmethod
    def save(cls, network=None, phases=[], filename=''):
        r"""
        Save all the pore and throat property data on the Network (and
        optionally on any Phases objects) to CSV files.

        Parameters
        ----------
        network : OpenPNM Network
            The Network containing the data to be stored

        phases : list of OpenPNM Phases (optional)
            The Phases whose data should be stored.

        filename : string or path object
            The name of the file to store the data

        Notes
        -----
        The data from all Geometry objects is added to the file automatically.

        """
        project, network, phases = cls._parse_args(network=network,
                                                   phases=phases)
        df = Pandas.to_dataframe(network=network, phases=phases, join=True)

        # Write to file
        if filename == '':
            filename = project.name
        fname = cls._parse_filename(filename=filename, ext='csv')
        df.to_csv(fname, index=False)

    @classmethod
    def load(cls, filename, project=None):
        r"""
        Opens a 'csv' file, reads in the data, and adds it to the **Network**

        Parameters
        ----------
        filename : string (optional)
            The name of the file containing the data to import.  The formatting
            of this file is outlined below.

        project : OpenPNM Project object
            A GenericNetwork is created and added to the specified Project.
            If no Project object is supplied then one will be created and
            returned.

        """
        net = {}

        fname = cls._parse_filename(filename)
        a = pd.read_table(filepath_or_buffer=fname,
                          sep=',',
                          skipinitialspace=True,
                          index_col=False,
                          true_values=['T', 't', 'True', 'true', 'TRUE'],
                          false_values=['F', 'f', 'False', 'false', 'FALSE'])

        dct = {}
        # First parse through all the items and clean-up`
        keys = sorted(list(a.keys()))
        for item in keys:
            # Merge arrays that have been split into multiple columns
            m = re.search(r'\[.\]', item)  # The dot '.' is a wildcard
            if m:  # m is None if pattern not found, otherwise merge cols
                pname = re.split(r'\[.\]', item)[0]  # Get base propname
                # Find all other keys with same base propname
                merge_keys = [k for k in a.keys() if k.startswith(pname)]
                # Rerieve and remove arrays with same base propname
                merge_cols = [a.pop(k) for k in merge_keys]
                # Merge arrays into multi-column array and store in DataFrame
                dct[pname] = sp.vstack(merge_cols).T
                # Remove key from list of keys
                [keys.pop(keys.index(k)) for k in keys if k.startswith(pname)]
            else:
                dct[item] = sp.array(a.pop(item))

        if project is None:
            project = Project(name=filename.split('.')[0])
        network = GenericNetwork(project=project)
        network = cls._update_network(network=network, net=net)
        return project
