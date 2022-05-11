from openpnm.io._pandas import Pandas
from openpnm.io import GenericIO
from openpnm.utils import Workspace

ws = Workspace()


class CSV(GenericIO):
    r"""
    Writes CSV (comma-separated-value files) containing pore and throat data

    """

    @classmethod
    def export_data(cls, network=None, phases=[], filename=''):
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
                                join=True, delim='.')
        # Write to file
        if filename == '':
            filename = project.name
        fname = cls._parse_filename(filename=filename, ext='csv')
        df.to_csv(fname, index=False)


def to_csv(network, phases=[], filename=''):
    CSV.export_data(network=network, phases=phases,
                    filename=filename)


to_csv.__doc__ = CSV.export_data.__doc__
