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
    def save(cls, simulation, filename=''):
        r"""
        Save all the pore and throat property data on the Network (and
        optionally on any Phases objects) to CSV files.

        Parameters
        ----------
        simulation : OpenPNM Simulation
            The Simulation containing the data to be stored

        filename : string
            The name of the file to store the data

        Notes
        -----
        The data from all Geometry objects is added to the file automatically.

        """
        dataframes = Pandas.get_data_frames(simulation=simulation)
        dfp = dataframes['pore.DataFrame']
        dft = dataframes['throat.DataFrame']
        b = dft.join(other=dfp, how='left')

        # Write to file
        if filename == '':
            filename = simulation.name
        with cls._write_file(filename=filename, ext='csv') as f:
            b.to_csv(f, index=False)

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
                              true_values=['T', 't', 'True', 'true',
                                           'TRUE'],
                              false_values=['F', 'f', 'False', 'false',
                                            'FALSE'])

        # Now parse through all the other items
        for item in a.keys():
            element = item.split('.')[0]
            prop = item.split('.', maxsplit=1)[1]
            data = sp.array(a[item].dropna())
            if type(data[0]) is str:
                N = sp.shape(data)[0]
                if '.' in data[0].split(' ')[0]:  # Decimal means float
                    dtype = float
                else:
                    dtype = int
                temp = sp.empty(sp.shape(data), dtype=object)
                for row in range(N):
                    temp[row] = sp.fromstring(data[row], sep=' ', dtype=dtype)
                data = sp.vstack(temp)
            else:
                dtype = type(data[0])
            net[element+'.'+prop] = data.astype(dtype)

        if simulation is None:
            simulation = Simulation(name=filename.split('.')[0])
        network = GenericNetwork(simulation=simulation)
        network = cls._update_network(network=network, net=net)
        return network
