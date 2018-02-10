import scipy as sp
import openpnm as op
from openpnm.core import logging
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

    @staticmethod
    def split_geometry(network):
        r"""
        This method accepts an OpenPNM Network object and removes all geometry
        related pore and throat properties, (basically all values other than
        ```'pore.coords'``` and ```throat.conns```), and places them on a
        GenericGeometry object.  Any labels on the Network are left intact.

        Parameters
        ----------
        network : OpenPNM Network Object
            The Network that possesses the geometrical values

        Returns
        -------
        geometry : OpenPNM Geometry Object
            The new GenericGeometry object that was created to contain the
            geometrical pore and throat properties.

        """
        geom = op.geometry.GenericGeometry(network=network,
                                           pores=network.Ps,
                                           throats=network.Ts)
        for item in network.props():
            if item not in ['pore.coords', 'throat.conns']:
                geom.update({item: network.pop(item)})
        return geom

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
    def _write_file(cls, filename, ext):
        ext = ext.replace('.', '').lower()
        filename = filename.rstrip('.'+ext)
        filename = filename+'.'+ext
        try:
            logger.warning(filename+' already exists, contents will be ' +
                           'overwritten')
            f = open(filename, mode='w')
        except:
            f = open(filename, mode='x')
        return f

    @classmethod
    def _read_file(cls, filename, ext, mode='r'):
        ext = ext.replace('.', '').lower()
        if not filename.endswith('.'+ext):
            filename = filename+'.'+ext
        f = open(filename, mode=mode)
        return f
