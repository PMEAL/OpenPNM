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
            interleave=True, flatten=True, categorize=False):
        r"""
        Returns a single dictionary object containing data from the given
        objects, with the keys organized differently depending on optional
        arguments.

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

        categorize : boolean (default is ``False``)
            If ``True`` the dictionary keys will be stored under a general
            level corresponding to their type (e.g. 'network/net_01/pore.all').
            If  ``interleave`` is ``True`` then only the only categories are
            *network* and *phase*, since *geometry* and *physics* data get
            stored under their respective *network* and *phase*.

        Returns
        -------
        A dictionary with all the data stored at the top level, but with the
        keys chosen to represent the hierarchical data structure.  The actual
        format of the keys depends on the arguments to the function.

        Notes
        -----
        The choice of '/' as a delimiter is chosen to work with the hdf5
        format.

        There is a handy package called *flatdict* that can be used to
        access this dictionary using normal dictionary syntax if so desired.

        """
        phases = cls._parse_phases(phases=phases)

        simulation = network.simulation
        d = PrintableDict()
        # This all still relies on automatic interleaving of data
        prefix = ''
        for key in network.keys(element=element):
            if categorize:
                prefix = 'network/'
            if interleave:
                d[prefix+network.name+'/'+key] = network[key]
            else:
                d[prefix+network.name+'/'+key] = network[key]
        for geo in simulation.geometries.values():
            for key in geo.keys(element=element):
                if interleave:
                    d[prefix+network.name+'/'+key] = network[key]
                else:
                    if flatten:
                        if categorize:
                            prefix = 'geometry/'
                        d[prefix+geo.name+'/'+key] = geo[key]
                    else:
                        d[prefix+network.name+'/'+geo.name+'/'+key] = geo[key]
        for phase in phases:
            for key in phase.keys(element=element):
                if categorize:
                    prefix = 'phase/'
                if interleave:
                    d[prefix+phase.name+'/'+key] = phase[key]
                else:
                    d[prefix+phase.name+'/'+key] = phase[key]
            for physics in simulation.find_physics(phase=phase):
                phys = simulation.physics[physics]
                for key in phys.keys(element=element):
                    if interleave:
                        d[prefix+phase.name+'/'+key] = phase[key]
                    else:
                        if flatten:
                            if categorize:
                                prefix = 'physics/'
                            d[prefix+phys.name+'/'+key] = phys[key]
                        else:
                            d[prefix+phase.name+'/'+phys.name+'/'+key] = phys[key]
        return d

    @classmethod
    def save(cls, network, phases=[], filename=''):
        r"""
        Saves data from the given objects into the specified file.

        Parameters
        ----------
        network : OpenPNM Network Object
            The network containing the desired data

        phases : list of OpenPNM Phase Objects (optional, default is none)
            A list of phase objects whose data are to be included

        Notes
        -----
        This function enforces a specific dictionary format, so that it
        can be consistently interpreted by the ``load`` function.  To get
        a dictionary with a different format use the ``get`` method, and then
        optionally save it to a file manually using the ``pickle`` standard
        library.

        This method only saves the data, not any of the pore-scale models or
        other attributes.  To save an actual OpenPNM Simulation use the
        ``Workspace`` object.

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
        Load data from the specified file into an OpenPNM simulation

        Parameters
        ----------
        filname : string
            The path to the file to be openned

        simulation : OpenPNM Simulation object
            A GenericNetwork is created and added to the specified Simulation.
            If no Simulation object is supplied then one will be created and
            returned.

        Notes
        -----
        This function is designed to open files creating using the ``save``
        function, which have a specific format.

        """
        with cls._read_file(filename=filename, ext='dct', mode='rb') as f:
            net = pickle.load(f)
        if simulation is None:
            simulation = Simulation(name=filename.split('.')[0])
        network = GenericNetwork(simulation=simulation)
        network = cls._update_network(network=network, net=net)
        return simulation
