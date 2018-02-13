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
    def to_hdf5(cls, network, phases=[], element=['pore', 'throat'],
                filename='', interleave=True, flatten=False,
                categorize_objects=False, categorize_data=False):
        r"""

        """
        dct = Dict.to_dict(network=network, phases=phases, element=element,
                           interleave=interleave, flatten=flatten,
                           categorize=categorize_objects)
        d = FlatDict(dct, delimiter='/')
        if filename == '':
            filename = network.simulation.name
        f = h5py.File(filename+".hdf5", "w")
        for item in d.keys():
            tempname = '_'.join(item.split('.'))
            arr = d[item]
            if 'U' in str(arr[0].dtype):
                pass
            else:
                if categorize_data:
                    temp = item.split('/')
                    if arr.dtype == bool:
                        temp = '/'.join(temp[:-1]) + '/label/' + temp[-1]
                        f.create_dataset(name=temp, shape=arr.shape,
                                         dtype=bool, data=arr)
                    else:
                        temp = '/'.join(temp[:-1]) + '/property/' + temp[-1]
                        f.create_dataset(name=temp, shape=arr.shape,
                                         dtype=arr.dtype, data=arr)
                else:
                    f.create_dataset(name='/'+tempname, shape=arr.shape,
                                     dtype=arr.dtype, data=arr)
        return f

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
        This method only saves the data, not any of the pore-scale models or
        other attributes.  To save an actual OpenPNM Simulation use the
        ``Workspace`` object.

        """
        simulation = network.simulation
        if filename == '':
            filename = simulation.name
        else:
            filename = filename.rsplit('.hdf5', 1)[0]
        f = cls.to_hdf5(network=network, phases=phases, interleave=True)
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
        def print_hdf5(f, p='', indent='―'):
            for item in f.keys():
                if hasattr(f[item], 'keys'):
                    p = print_hdf5(f[item], p=p, indent=indent + indent[0])
                elif indent[-1] != ' ':
                    indent = indent + '| '
                p = indent + item + '\n' + p
            return(p)
        p = print_hdf5(f)
        print(p)
