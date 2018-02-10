import pickle
import scipy as sp
from openpnm.core import logging
from openpnm.network import GenericNetwork
from openpnm.utils import PrintableDict
from openpnm.io import GenericIO
logger = logging.getLogger(__name__)


class Dict(GenericIO):
    r"""

    """

    @classmethod
    def get(cls, simulation, phases=[]):
        r"""

        """
        if type(phases) is not list:  # Ensure it's a list
            phases = [phases]

        network = simulation.network
        d = PrintableDict()
        for key in network.keys():
            d[network.name+'.'+key] = network[key]
        for geo in simulation.geometries.values():
            for key in geo.keys():
                if network.name+'.'+key not in d.keys():
                    # This relies on automatic interleaving of data
                    d[network.name+'.'+key] = network[key]
        for phase in phases:
            for key in phase.keys():
                d[phase.name+'.'+key] = phase[key]
            for phys in simulation.find_physics(phase=phase):
                physics = simulation.physics[phys]
                for key in physics.keys():
                    if phase.name+'.'+key not in d.keys():
                        # This relies on automatic interleaving of data
                        d[phase.name+'.'+key] = phase[key]
        return d

    @classmethod
    def save(cls, simulation, filename=''):
        r"""

        """
        if filename == '':
            filename = simulation.name
        else:
            filename = filename.rsplit('.dct', 1)[0]
        d = cls.get(simulation=simulation, phases=[])
        for item in list(d.keys()):
            new_name = item.split('.')
            d[new_name[1] + '.' + new_name[2]] = d.pop(item)
        pickle.dump(d, open(filename + '.dct', 'wb'))

    @classmethod
    def load(cls, filename):
        r"""

        """
        with cls._read_file(filename=filename, ext='dct', mode='rb') as f:
            net = pickle.load(f)

        network = GenericNetwork()
        network = cls._update_network(network=network, net=net)
        return network.simulation
