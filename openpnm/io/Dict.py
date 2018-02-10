import pickle
from openpnm.core import logging, Simulation
from openpnm.network import GenericNetwork
from openpnm.utils import FlatDict, PrintableDict
from openpnm.io import GenericIO
logger = logging.getLogger(__name__)


class Dict(GenericIO):
    r"""

    """

    @classmethod
    def get(cls, network, phases=[], element=['pore', 'throat'],
            interleave=True, flatten=True):
        r"""
        Returns a single dictionary object containing data from the given
        objects.  The dictionary keys are organized differently depending on
        optional arguments.

        Parameters
        ----------
        network : OpenPNM Network Object
            The network containing the desired data

        phases : list of OpenPNM Phase Objects (optional, default is none)
            A list of phase objects whose data are to be included

        element : string or list of strings
            An indication of whether 'pore' and/or 'throat' data are desired.
            The default is both.

        interleave : boolean (default is ``True``)
            When ``True`` (default) the data from all Geometry objects (and
            Physics objects if ``phases`` are given) is interleaved into
            a single array and stored as a network property (or Phase
            property for Physics data). When ``False``, the data for each
            object are stored under their own dictionary key, the structuring
            of which depends on the value of the ``flatten`` argument.

        flatten : boolean (default is ``True``)
            When ``True``, all objects are accessible from the top level
            of the dictionary.  When ``False`` objects are nested under their
            parent object.  If ``interleave`` is ``True`` this argument is
            ignored.

        Returns
        -------
        This returns a special dictionary type called a FlatDict (from the
        *flatdict* package), that is designed to access hierarchical (i.e.
        nested) dictionaries .  A FlatDict can use normal dict syntax such as:

        ``>>> array = flatdict['top_dict']['inner_dict']['value']``

        but can also use a single key syntax such as:

        ``>>> array = flatdict['top_dict/inner_dict/value']``

        Importantly, the ``keys`` attribute returns the latter values so one
        can directly visualize the hierarchical structure.

        Notes
        -----
        The choice of '/' as a delimiter is chosen to work with the hdf5 format

        """
        if type(phases) is not list:  # Ensure it's a list
            phases = [phases]

        simulation = network.simulation
        d = FlatDict(delimiter='/')
        # This all still relies on automatic interleaving of data
        for key in network.keys(element=element):
            if interleave:
                d[network.name+'/'+key] = network[key]
            else:
                d[network.name+'/'+key] = network[key]
        for geo in simulation.geometries.values():
            for key in geo.keys(element=element):
                if interleave:
                    d[network.name+'/'+key] = network[key]
                else:
                    if flatten:
                        d[geo.name+'/'+key] = geo[key]
                    else:
                        d[network.name+'/'+geo.name+'/'+key] = geo[key]
        for phase in phases:
            for key in phase.keys(element=element):
                if interleave:
                    d[phase.name+'/'+key] = phase[key]
                else:
                    d[phase.name+'/'+key] = phase[key]
            for physics in simulation.find_physics(phase=phase):
                phys = simulation.physics[physics]
                for key in phys.keys(element=element):
                    if interleave:
                        d[phase.name+'/'+key] = phase[key]
                    else:
                        if flatten:
                            d[phys.name+'/'+key] = phys[key]
                        else:
                            d[phase.name+'/'+phys.name+'/'+key] = phys[key]
        return d

    @classmethod
    def save(cls, network, phases=[], filename=''):
        r"""

        """
        simulation = network.simulation
        if filename == '':
            filename = simulation.name
        else:
            filename = filename.rsplit('.dct', 1)[0]
        d = cls.get(simulation=simulation, phases=phases, interleave=True)
        for item in list(d.keys()):
            new_name = item.split('.')
            d[new_name[1] + '.' + new_name[2]] = d.pop(item)
        pickle.dump(d, open(filename + '.dct', 'wb'))

    @classmethod
    def load(cls, filename, simulation=None):
        r"""

        """
        with cls._read_file(filename=filename, ext='dct', mode='rb') as f:
            net = pickle.load(f)
        if simulation is None:
            simulation = Simulation(name=filename.split('.')[0])
        network = GenericNetwork(simulation=simulation)
        network = cls._update_network(network=network, net=net)
        return simulation
