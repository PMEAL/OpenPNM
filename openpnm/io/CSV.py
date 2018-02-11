import re
import scipy as sp
import pandas as pd
from openpnm.core import logging, Simulation
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
    def save(cls, network, phases=[], filename=''):
        r"""
        Save all the pore and throat property data on the Network (and
        optionally on any Phases objects) to CSV files.

        Parameters
        ----------
        network : OpenPNM Network
            The Network containing the data to be stored

        phases : list of OpenPNM Phases (optional)
            The Phases whose data should be stored.

        filename : string
            The name of the file to store the data

        Notes
        -----
        The data from all Geometry objects is added to the file automatically.

        """
        phases = cls._parse_phases(phases=phases)
        simulation = network.simulation

        df = Pandas.get(network=network, phases=phases, join=True)

        # Write to file
        if filename == '':
            filename = simulation.name
        with cls._write_file(filename=filename, ext='csv') as f:
            df.to_csv(f, index=False)

    @classmethod
    def load(cls, filename, simulation=None):
        r"""
        Opens a 'csv' file, reads in the data, and adds it to the **Network**

        Parameters
        ----------
        filename : string (optional)
            The name of the file containing the data to import.  The formatting
            of this file is outlined below.

        simulation : OpenPNM Simulation object
            A GenericNetwork is created and added to the specified Simulation.
            If no Simulation object is supplied then one will be created and
            returned.

        """
        net = {}

        with cls._read_file(filename=filename, ext='csv') as f:
            a = pd.read_table(filepath_or_buffer=f,
                              sep=',',
                              skipinitialspace=True,
                              index_col=False,
                              true_values=['T', 't', 'True', 'true', 'TRUE'],
                              false_values=['F', 'f', 'False', 'false',
                                            'FALSE'])

        net = {}
        # Now parse through all the items and clean-up
        pat = r'\[.\]'  # This pattern is used to find columns with indices
        keys = list(a.keys())  # Use a list so it can be mutated inside loop
        for item in keys:
            # Deal with arrays that have been split into multiple columns
            m = re.search(pat, item)
            if m:  # m is None if pattern not found, otherwise merge occurences
                pname = re.split(pat, item)[0]  # Get base propname
                # Find all other keys with same base propname
                all_keys = sorted([k for k in keys if k.startswith(pname)])
                # Rerieve and remove arrays with same base propname
                all_arrays = [a.pop(k) for k in all_keys]
                # Remove other keys with same base propname
                keys = [keys.pop(keys.index(k)) for k in all_keys]
                # Merge arrays into multi-column array and store in DataFrame
                net[pname] = sp.vstack(all_arrays).T
            else:
                net[item] = sp.array(a[item].dropna())

        if simulation is None:
            simulation = Simulation(name=filename.split('.')[0])
        network = GenericNetwork(simulation=simulation)
        network = cls._update_network(network=network, net=net)
        return simulation
