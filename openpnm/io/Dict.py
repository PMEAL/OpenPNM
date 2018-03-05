import pickle
from openpnm.core import logging, Project, Workspace
from openpnm.network import GenericNetwork
from openpnm.utils import NestedDict, FlatDict
from openpnm.io import GenericIO
logger = logging.getLogger(__name__)
ws = Workspace()


class Dict(GenericIO):
    r"""
    This is the most important class in the ``io`` module.  It is used to
    generate hierarchical ``dicts`` in a given format for use in many of the
    the other classes (e.g. ``Pandas`` and ``HDF5``, which are subsequently
    used in ``CSV`` and ``XDMF``, and so forth).

    """

    @classmethod
    def from_dict(cls, dct, project=None):
        r"""
        """
        if project is None:
            project = ws.new_project()
        dct = NestedDict(dct, delimiter=' | ')
        obj_types = ['network', 'geometry', 'phase', 'physics', 'algorithm']
        for item in dct.keys():
            if item in obj_types:
                for name in dct[item].keys():
                    try:
                        obj = project[name]
                    except KeyError:
                        obj = project._new_object(objtype=item, name=name)
                    obj.update(dct[item][name])
            else:
                name = item
                try:
                    obj = project[name]
                except KeyError:
                    obj = project._new_object(objtype='base', name=name)
                obj.update(dct[name])
        return project

    @classmethod
    def to_dict(cls, network=None, phases=[], element=['pore', 'throat'],
                interleave=True, flatten=True, categorize_by=[]):
        r"""
        Returns a single dictionary object containing data from the given
        OpenPNM objects, with the keys organized differently depending on
        optional arguments.

        Parameters
        ----------
        network : OpenPNM Network Object (optional)
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

        categorize_by : string or list of strings
            Indicates how the dictionaries should be organized.  The list can
            contain any, all or none of the following strings:

            **'object'** : If specified the dictionary keys will be stored
            under a general level corresponding to their type (e.g.
            'network/net_01/pore.all'). If  ``interleave`` is ``True`` then
            only the only categories are *network* and *phase*, since
            *geometry* and *physics* data get stored under their respective
            *network* and *phase*.

            **'data'** : If specified the data arrays are additionally
            categorized by ``label`` and ``property`` to separate *boolean*
            from *numeric* data.

            **'element'** : If specified the data arrays are
            additionally categorized by ``pore`` and ``throat``, meaning
            that the propnames are no longer prepended by a 'pore.' or
            'throat.'

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
        project, network, phases = cls._parse_args(network=network,
                                                   phases=phases)

        d = NestedDict()
        delim = '/'

        def build_path(obj, key):
            propname = delim + key
            prefix = 'root'
            datatype = ''
            arr = obj[key]
            if 'object' in categorize_by:
                prefix = obj._isa()
            if 'element' in categorize_by:
                propname = delim + key.replace('.', delim)
            if 'data' in categorize_by:
                if arr.dtype == bool:
                    datatype = delim + 'labels'
                else:
                    datatype = delim + 'properties'
            path = prefix + delim + obj.name + datatype + propname
            return path

        for net in network:
            for key in net.keys(element=element, mode='all'):
                path = build_path(obj=net, key=key,)
                d[path] = net[key]

            for geo in project.geometries().values():
                for key in geo.keys(element=element, mode='all'):
                    if interleave:
                        path = build_path(obj=net, key=key)
                        d[path] = net[key]
                    else:
                        path = build_path(obj=geo, key=key)
                        if flatten:
                            d[path] = geo[key]
                        elif 'object' not in categorize_by:
                            path = path.split(delim)
                            path.insert(1, net.name)
                            path = delim.join(path)
                        d[path] = geo[key]

        for phase in phases:
            for key in phase.keys(element=element, mode='all'):
                path = build_path(obj=phase, key=key)
                d[path] = phase[key]

            for phys in project.find_physics(phase=phase):
                for key in phys.keys(element=element, mode='all'):
                    if interleave:
                        path = build_path(obj=phase, key=key)
                        d[path] = phase[key]
                    else:
                        path = build_path(obj=phys, key=key)
                        if flatten:
                            d[path] = phys[key]
                        elif 'object' not in categorize_by:
                            path = path.split(delim)
                            path.insert(1, phase.name)
                            path = delim.join(path)
                        d[path] = phys[key]

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
        other attributes.  To save an actual OpenPNM Project use the
        ``Workspace`` object.

        """
        project = network.project
        if filename == '':
            filename = project.name
        else:
            filename = filename.rsplit('.dct', 1)[0]
        d = cls.to_dict(project=project, phases=phases,
                        interleave=True, categorize=False)
        for item in list(d.keys()):
            new_name = item.split('.')
            d[new_name[1] + '.' + new_name[2]] = d.pop(item)
        pickle.dump(d, open(filename + '.dct', 'wb'))

    @classmethod
    def load(cls, filename, project=None):
        r"""
        Load data from the specified file into an OpenPNM project

        Parameters
        ----------
        filname : string
            The path to the file to be openned

        project : OpenPNM Project object
            A GenericNetwork is created and added to the specified Project.
            If no Project object is supplied then one will be created and
            returned.

        Notes
        -----
        This function is designed to open files creating using the ``save``
        function, which have a specific format.

        """
        with cls._read_file(filename=filename, ext='dct', mode='rb') as f:
            net = pickle.load(f)
        if project is None:
            project = Project(name=filename.split('.')[0])
        network = GenericNetwork(project=project)
        network = cls._update_network(network=network, net=net)
        return project
