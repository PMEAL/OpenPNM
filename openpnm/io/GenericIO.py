import scipy as sp
import numpy as np
import openpnm as op
from pathlib import Path
from openpnm.core import logging
from openpnm.utils import flat_list
logger = logging.getLogger(__name__)


class GenericIO():

    @classmethod
    def save(cls):
        raise NotImplementedError("The \'save\' method for this class " +
                                  "does not exist yet")

    @classmethod
    def load(cls):
        raise NotImplementedError("The \'load\' method for this class " +
                                  "does not exist yet")

    @classmethod
    def _convert_data(cls, project):
        # Convert arrays of 1's and/or 0's to booleans
        for obj in project:
            for item in obj.keys():
                ra = obj[item]
                if (ra.dtype == int) and (ra.max() <= 1) and (ra.min() >= 0):
                    obj.update({item: ra.astype(bool)})
        return project

    @classmethod
    def _update_network(cls, network, net):
        # Infer Np and Nt from length of given prop arrays in file
        for el in ['pore', 'throat']:
            N = [sp.shape(net[i])[0] for i in net.keys() if i.startswith(el)]
            if N:
                N = sp.array(N)
                if sp.all(N == N[0]):
                    if (network._count(el) == N[0]) \
                            or (network._count(el) == 0):
                        network.update({el+'.all': sp.ones((N[0],),
                                                           dtype=bool)})
                        net.pop(el+'.all', None)
                    else:
                        raise Exception('Length of '+el+' data in file' +
                                        ' does not match network')
                else:
                    raise Exception(el+' data in file have inconsistent' +
                                    ' lengths')
        # Add data on dummy net to actual network
        for item in net.keys():
            # Try to infer array types and change if necessary
            # Chcek for booleans disguised and 1's and 0's
            if (net[item].dtype is int) and (net[item].max() == 1):
                net[item] = net[item].astype(bool)
            # Write data to network object
            if item in network:
                logger.warning('\''+item+'\' already present...overwriting')
            network.update({item: net[item]})

        network._gen_ids()

        return network

    @classmethod
    def _parse_filename(cls, filename, path='', ext=''):
        p = Path(path)
        p = p.resolve()
        p = p.joinpath(filename)
        # If extension not part of filename
        if p.suffix == '':
            try:
                p = p.with_suffix(ext)
            except ValueError:
                p = p.with_suffix('.'+ext)
        return p

    @classmethod
    def _parse_args(cls, network, phases):
        # Convert network to a list, even if empty
        if network is None:
            network = []
        else:
            network = flat_list(network)
        # Ensure phases is a list, even if empty
        phases = flat_list(phases)
        # Get handle to project object
        if len(network) == 0:
            if len(phases) == 0:
                raise Exception('Must specify one of network or phase')
            else:
                project = phases[0].project
        else:
            project = network[0].project
        return (project, network, phases)
