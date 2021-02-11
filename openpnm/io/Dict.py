import pickle
from flatdict import FlatDict
from openpnm.utils import NestedDict, sanitize_dict, Workspace
from openpnm.utils import logging
from openpnm.io import GenericIO
logger = logging.getLogger(__name__)
ws = Workspace()


class Dict(GenericIO):
    r"""
    Generates hierarchical ``dicts`` with a high degree of control over the
    structure.

    This is the most important class in the ``io`` module, since many other
    classes use this to manipulate and format the data structures.

    Also, it is possible to use Python's ``pickle`` module to save ``dicts``
    to file.

    """

    @classmethod
    def from_dict(cls, dct, project=None, delim=' | '):
        r"""
        This method converts a correctly formatted dictionary into OpenPNM
        objects, and returns a handle to the *project* containing them.

        Parameters
        ----------
        dct : dictionary
            The Python dictionary containing the data.  The nesting and
            labeling of the dictionary is used to create the appropriate
            OpenPNM objects.

        project : OpenPNM Project Object
            The project with which the created objects should be associated.
            If not supplied, one will be created.

        Returns
        -------
        project : list
            An OpenPNM Project containing the objects created to store the
            given data.

        Notes
        -----
        The requirement of a *correctly formed* dictionary is rather strict,
        and essentially means a dictionary produced by the ``to_dict`` method
        of this class.

        """
        if project is None:
            project = ws.new_project()

        # Uncategorize pore/throat and labels/properties, if present
        fd = FlatDict(dct, delimiter=delim)
        # If . is the delimiter, replace with | otherwise things break
        if delim == '.':
            delim = ' | '
            for key in list(fd.keys()):
                new_key = key.replace('.', delim)
                fd[new_key] = fd.pop(key)
        d = FlatDict(delimiter=delim)
        for key in list(fd.keys()):
            new_key = key.replace('pore' + delim, 'pore.')
            new_key = new_key.replace('throat' + delim, 'throat.')
            new_key = new_key.replace('labels' + delim, '')
            new_key = new_key.replace('properties' + delim, '')
            d[new_key] = fd.pop(key)

        # Plase data into correctly categorized dicts, for later handling
        objs = {'network': NestedDict(),
                'geometry': NestedDict(),
                'physics': NestedDict(),
                'phase': NestedDict(),
                'algorithm': NestedDict(),
                'base': NestedDict()}
        for item in d.keys():
            path = item.split(delim)
            if len(path) > 2:
                if path[-3] in objs.keys():
                    # Item is categorized by type, so note it
                    objs[path[-3]][path[-2]][path[-1]] = d[item]
                else:
                    # Item is nested, not categorized; make it a base
                    objs['base'][path[-2]][path[-1]] = d[item]
            else:
                # If not categorized by type, make it a base
                objs['base'][path[-2]][path[-1]] = d[item]

        # Convert to OpenPNM Objects, attempting to infer type
        for objtype in objs.keys():
            for name in objs[objtype].keys():
                # Create empty object, using dummy name to avoid error
                obj = project._new_object(objtype=objtype, name='')
                # Overwrite name
                obj._set_name(name=name, validate=False)
                # Update new object with data from dict
                obj.update(objs[objtype][name])

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
        delim = ' | '
        d = NestedDict(delimiter=delim)

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
                path = build_path(obj=net, key=key)
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
                        elif 'object' in categorize_by:
                            path = path.split(delim)
                            path.insert(0, 'network')
                            path.insert(1, net.name)
                            path = delim.join(path)
                        else:
                            path = path.split(delim)
                            path.insert(1, net.name)
                            path = delim.join(path)
                        d[path] = geo[key]

        for phase in phases:
            for key in phase.keys(element=element, mode='all'):
                path = build_path(obj=phase, key=key)
                d[path] = phase[key]

            for phys in project.find_physics(phase=phase):
                if phys:
                    for key in phys.keys(element=element, mode='all'):
                        if interleave:
                            path = build_path(obj=phase, key=key)
                            d[path] = phase[key]
                        else:
                            path = build_path(obj=phys, key=key)
                            if flatten:
                                d[path] = phys[key]
                            elif 'object' in categorize_by:
                                path = path.split(delim)
                                path.insert(0, 'phase')
                                path.insert(1, phase.name)
                                path = delim.join(path)
                            else:
                                path = path.split(delim)
                                path.insert(1, phase.name)
                                path = delim.join(path)
                            d[path] = phys[key]

        if 'root' in d.keys():
            d = d['root']
        if 'project' in categorize_by:
            new_d = NestedDict()
            new_d[project.name] = d
            d = new_d
        return d

    @classmethod
    def save(cls, *args, **kwargs):
        r"""
        This method is being deprecated.  Use ``export_data`` instead.
        """
        cls.export_data(*args, **kwargs)

    @classmethod
    def export_data(cls, dct, filename):
        r"""
        Saves data from the given dictionary into the specified file.

        Parameters
        ----------
        dct : dictionary
            A dictionary to save to file, presumably obtained from the
            ``to_dict`` method of this class.

        filename : string or path object
            The filename to store the dictionary.

        Notes
        -----
        This method uses the pickle module to save the dictionary.

        """
        fname = cls._parse_filename(filename=filename, ext='dct')
        dct = sanitize_dict(dct)
        with open(fname, 'wb') as f:
            pickle.dump(dct, f)

    @classmethod
    def load(cls, *args, **kwargs):
        r"""
        This method is being deprecated.  Use ``import_data`` instead.
        """
        return cls.import_data(*args, **kwargs)

    @classmethod
    def import_data(cls, filename):
        r"""
        Load data from the specified pickle file into a Python dictionary

        Parameters
        ----------
        filename : string
            The path to the file to be opened

        Notes
        -----
        This returns a Python dictionary which can be converted into OpenPNM
        objects using the ``from_dict`` method of this class.

        """
        fname = cls._parse_filename(filename)
        with open(fname, 'rb') as f:
            dct = pickle.load(f)
        return dct
