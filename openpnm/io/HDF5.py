import h5py
from openpnm.core import logging, Simulation
from openpnm.io import Dict
from openpnm.network import GenericNetwork
from openpnm.utils import FlatDict
from openpnm.io import GenericIO
logger = logging.getLogger(__name__)


class HDF5(GenericIO):
    r"""

    """

    @classmethod
    def to_hdf5(cls, network=None, phases=[], element=['pore', 'throat'],
                filename='', interleave=True, flatten=False, categorize_by=[]):
        r"""
        Creates an HDF5 file containing data from the specified objects,
        and categorized according to the given arguments.

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

        categorize_by : string or list of strings
            Indicates how the dictionaries should be organized.  The list can
            contain any, all or none of the following strings:

            **'objects'** : If specified the dictionary keys will be stored
            under a general level corresponding to their type (e.g.
            'network/net_01/pore.all'). If  ``interleave`` is ``True`` then
            only the only categories are *network* and *phase*, since
            *geometry* and *physics* data get stored under their respective
            *network* and *phase*.

            **'data'** : If specified the data arrays are additionally
            categorized by ``label`` and ``property`` to separate *boolean*
            from *numeric* data.

            **'categorize_elements'** : If specified the data arrays are
            additionally categorized by ``pore`` and ``throat``, meaning
            that the propnames are no longer prepended by a 'pore.' or
            'throat.'

        """
        simulation, network, phases = cls._parse_args(network=network,
                                                      phases=phases)
        dct = Dict.to_dict(network=network, phases=phases, element=element,
                           interleave=interleave, flatten=flatten,
                           categorize_by=categorize_by)
        d = FlatDict(dct, delimiter='/')
        if filename == '':
            filename = simulation.name
        f = h5py.File(filename+".hdf", "w")
        for item in d.keys():
            tempname = '_'.join(item.split('.'))
            arr = d[item]
            if 'U' in str(arr[0].dtype):
                pass
            else:
                f.create_dataset(name='/'+tempname, shape=arr.shape,
                                 dtype=arr.dtype, data=arr)
        return f

    @classmethod
    def save(cls, network=None, phases=[], filename='', **kwargs):
        r"""
        Saves data from the given objects into the specified file.

        Parameters
        ----------
        network : OpenPNM Network Object
            The network containing the desired data

        phases : list of OpenPNM Phase Objects (optional, default is none)
            A list of phase objects whose data are to be included

        filename : string

        **kwargs : key-word arguments
            These additional arguments are passed on to the ``to_hdf5``
            function to controls the hierarchy of the data.  Refer to that
            method's docstring for more info.

        Notes
        -----
        This method only saves the data, not any of the pore-scale models or
        other attributes.  To save an actual OpenPNM Simulation use the
        ``Workspace`` object.

        """
        simulation, network, phases = cls._parse_args(network=network,
                                                      phases=phases)

        if filename == '':
            filename = simulation.name
        else:
            filename = filename.rsplit('.hdf', 1)[0]
        f = cls.to_hdf5(network=network, phases=phases, interleave=True,
                        **kwargs)
        f.close()

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
        raise NotImplementedError()

    def print_hierarchy(f):
        def print_level(f, p='', indent='-'):
            for item in f.keys():
                if hasattr(f[item], 'keys'):
                    p = print_level(f[item], p=p, indent=indent + indent[0])
                elif indent[-1] != ' ':
                    indent = indent + ''
                p = indent + item + '\n' + p
            return(p)
        p = print_level(f)
        print(p)

    def print_flattened(f):
        f.visit(print)
