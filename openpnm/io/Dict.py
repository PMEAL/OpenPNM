import pickle
from openpnm.core import logging, Simulation, Base
from openpnm.network import GenericNetwork
from openpnm.utils import NestedDict, FlatDict
from openpnm.io import GenericIO
logger = logging.getLogger(__name__)


class Dict(GenericIO):
    r"""
    This is the most important class in the ``io`` module.  It is used to
    generate hierarchical ``dicts`` in a given format for use in many of the
    the other classes (e.g. ``Pandas`` and ``HDF5``, which are subsequently
    used in ``CSV`` and ``XDMF``, and so forth).

    """

    @classmethod
    def from_dict(cls, dct, simulation=None):

        # Now parse through dict and put values on correct objects
        dct = FlatDict(dct, delimiter='/')
        sim = NestedDict()
        for item in dct.keys():
            level = item.split('/')
            sim[level[-2]][level[-1]] = dct[item]

        if simulation is None:
            simulation = Simulation()

        for item in sim.keys():
            obj = Base(simulation=simulation)
            obj.update(sim[item])
            obj._name = item

        return simulation

    @classmethod
    def to_dict(cls, network, phases=[], element=['pore', 'throat'],
                interleave=True, flatten=True, categorize_objects=False,
                categorize_data=False):
        r"""
        Returns a single dictionary object containing data from the given
        OpenPNM objects, with the keys organized differently depending on
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

        categorize_objects : boolean (default is ``False``)
            If ``True`` the dictionary keys will be stored under a general
            level corresponding to their type (e.g. 'network/net_01/pore.all').
            If  ``interleave`` is ``True`` then only the only categories are
            *network* and *phase*, since *geometry* and *physics* data get
            stored under their respective *network* and *phase*.

        categorize_data : boolean (default is ``False)
            If ``True`` the data arrays are additionally categorized by
            ``label`` and ``property`` to separate *boolean* from *numeric*
            data.

        categorize_elements : boolean (default is ``False)
            If ``True`` the data arrays are additionally categorized by
            ``pore`` and ``throat``.

        Returns
        -------
        A dictionary with the data stored in a hierarchical data structure, the
        actual format of which depends on the arguments to the function.

        Notes
        -----
        There is a handy package called *flatdict* that can be used to
        access this dictionary using a single key such that:

        ``d[level_1][level_2] == d[level_1/level_2]``

        Importantly, converting to a *flatdict* allows it be converted to an
        *HDF5* file directly, since the hierarchy is dictated by the placement
        of '/' characters.
        """
        phases = cls._parse_phases(phases=phases)

        simulation = network.simulation
        d = NestedDict()
        # This all still relies on automatic interleaving of data
        prefix = 'root'
        for key in network.keys(element=element):
            if categorize_objects:
                prefix = 'network'
            if interleave:
                d[prefix][network.name][key] = network[key]
            else:
                d[prefix][network.name][key] = network[key]
        for geo in simulation.geometries().values():
            for key in geo.keys(element=element):
                if interleave:
                    d[prefix][network.name][key] = network[key]
                else:
                    if flatten:
                        if categorize_objects:
                            prefix = 'geometry'
                        d[prefix][geo.name][key] = geo[key]
                    else:
                        d[prefix][network.name][geo.name][key] = geo[key]
        for phase in phases:
            for key in phase.keys(element=element):
                if categorize_objects:
                    prefix = 'phase'
                if interleave:
                    d[prefix][phase.name][key] = phase[key]
                else:
                    d[prefix][phase.name][key] = phase[key]
            for physics in simulation.find_physics(phase=phase):
                phys = simulation.physics()[physics]
                for key in phys.keys(element=element):
                    if interleave:
                        d[prefix][phase.name][key] = phase[key]
                    else:
                        if flatten:
                            if categorize_objects:
                                prefix = 'physics'
                            d[prefix][phys.name][key] = phys[key]
                        else:
                            d[prefix][phase.name][phys.name][key] = phys[key]
        if 'root' in d.keys():
            d = d['root']
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
        d = cls.to_dict(simulation=simulation, phases=phases,
                        interleave=True, categorize=False)
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
