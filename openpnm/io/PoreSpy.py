import pickle as pk
from openpnm.utils import logging
from openpnm.io import GenericIO
from openpnm.network import GenericNetwork
from openpnm.geometry import Imported
logger = logging.getLogger(__name__)


class PoreSpy(GenericIO):
    r"""
    """

    @classmethod
    def load(cls, filename, project=None, settings={}):
        r"""
        """
        # Parse the link1 file
        filename = cls._parse_filename(filename=filename, ext='dict')
        with open(filename, mode='rb') as f:
            net = pk.load(f)

        network = GenericNetwork(project=project)
        network = cls._update_network(network=network, net=net)
        # The following list of props comes from porespy but are not good
        for item in ['pore.diameter',
                     'pore.area',
                     'throat.conduit_lengths.pore1',
                     'throat.conduit_lengths.pore2',
                     'throat.conduit_lengths.throat',
                     'throat.endpoints.head',
                     'throat.endpoints.tail',
                     'throat.diameter',
                     'throat.length',
                     'throat.volume']:
            network.pop(item, None)

        # Deal with a few other PoreSpy specific issues
        network['pore.region_volume'] = network.pop('pore.volume')

        Imported(network=network, settings=settings)

        return network.project
