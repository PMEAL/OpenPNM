import pickle as pk
from openpnm.io import GenericIO
from openpnm.network import GenericNetwork
from openpnm.geometry import Imported


class PoreSpy(GenericIO):
    r"""
    """

    @classmethod
    def import_data(cls, filename, project=None, settings=None):
        r"""
        Load a network extracted using the PoreSpy package

        Parameters
        ----------
        filename : str or dict
            Can either be a filename point to a pickled dictionary, or an
            actual dictionary.  The second option lets users avoid the
            step of saving the dictionary to a file
        project : Project
            If given, the loaded network and geometry will be added to this
            project, otherwise a new one will be created.

        """
        # Parse the filename
        if isinstance(filename, dict):
            net = filename
        else:
            filename = cls._parse_filename(filename=filename)
            with open(filename, mode='rb') as f:
                net = pk.load(f)

        network = GenericNetwork(project=project)
        network = cls._update_network(network=network, net=net)
        Imported(network=network, settings=settings)

        return network.project


def from_porespy(filename, project=None, settings=None):
    project = PoreSpy.import_data(filename=filename, project=project, settings=settings)
    return project


from_porespy.__doc__ = PoreSpy.import_data.__doc__
