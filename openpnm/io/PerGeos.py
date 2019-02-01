import scipy as sp
from openpnm.utils import logging
from openpnm.io import GenericIO
from openpnm.network import GenericNetwork
logger = logging.getLogger(__name__)


class PerGeos(GenericIO):
    r"""

    """

    @classmethod
    def save(cls, network=None, phases=[], filename=''):
        r"""

        """
        project, network, phases = cls._parse_args(network=network,
                                                   phases=phases)

        s = ["# Avizo 3D ASCII 3.0\n\n"]
        s.append("define VERTEX " + str(network.Np) + '\n')
        s.append("define EDGE " + str(network.Nt) + '\n')
        s.append("define VERTEX " + str(2*network.Nt) + '\n')
        s.append("Parameters {\n\tContentType \"HxPoreNetworkModel\"\n}\n\n")

        types = {'b': 'int', 'i': 'int', 'f': 'float'}
        typemap = {}
        namemap = {}
        shapemap = {}
        propmap = {}
        i = 1
        for item in network.keys():
            typemap[item] = types[str(network[item].dtype)[0]]
            ncols = int(network[item].size/network[item].shape[0])
            if ncols > 1:
                shapemap[item] = '[' + str(ncols) + ']'
            else:
                shapemap[item] = ''
            if item.startswith('pore'):
                element = 'pore', 'VERTEX'
            if item.startswith('throat'):
                element = 'throat', 'EDGE'
            n = item.replace(element[0] + '.', '').replace('.', '_').split('_')
            n = ''.join([i[0].upper()+i[1:] for i in n if len(i)])
            namemap[item] = n
            temp = element[1] + " { " + typemap[item] + shapemap[item] + " " +\
                   namemap[item] + " } @" + str(i) + '\n'
            s.append(temp)
            i += 1
            propmap[item] = str(i)

        s.append("# Data section follows")
        for item in network.keys():
            data = network[item]
            s.append('\n\n@' + propmap[item] + '\n')
            if shapemap[item] == '':
                data = sp.atleast_2d(data).T
            if typemap[item] == 'float':
                formatter = {'float_kind': lambda x: "%.15E" % x}
            else:
                formatter = None
            if data.dtype == 'bool':
                data = data.astype(int)
            d = sp.array2string(data, formatter=formatter)
            s.append(d.replace('[', '').replace(']', '').replace('\n ', '\n'))

        # Write to file
        if filename == '':
            filename = project.name
        fname = cls._parse_filename(filename=filename, ext='am')
        with open(fname, 'w') as f:
            f.write(''.join(s))

    @classmethod
    def load(cls, filename, network=None):
        r"""


        """
        net = {}

        # ---------------------------------------------------------------------
        # Parse the link1 file
        filename = cls._parse_filename(filename=filename, ext='am')
        with open(filename, mode='r') as f:
            Np = None
            Nt = None
            while (Np is None) or (Nt is None):
                s = f.readline()[:-1].split(' ')
                if s[0] == 'define':
                    if s[1] == 'VERTEX':
                        Np = int(s[2])
                    if s[1] == 'EDGE':
                        Nt = int(s[2])

            net = {}
            propmap = {}
            typemap = {}
            shapemap = {}
            while True:
                s = f.readline()[:-1].split(' ')
                if s[0] == 'VERTEX':
                    dshape = [Np]
                    if s[2].endswith(']'):
                        ncols = int(s[2].split('[', 1)[1].split(']')[0])
                        dshape.append(ncols)
                    dtype = s[2].split('[')[0]
                    temp = sp.zeros(dshape, dtype=dtype)
                    net['pore.'+s[3]] = temp
                    key = int(s[-1].replace('@', ''))
                    propmap[key] = 'pore.'+s[3]
                    typemap[key] = dtype
                    shapemap[key] = dshape
                elif s[0] == 'EDGE':
                    dshape = [Nt]
                    if s[2].endswith(']'):
                        ncols = int(s[2].split('[', 1)[1].split(']')[0])
                        dshape.append(ncols)
                    dtype = s[2].split('[')[0]
                    temp = sp.zeros(dshape, dtype=dtype)
                    net['throat.'+s[3]] = temp
                    key = int(s[-1].replace('@', ''))
                    propmap[key] = 'throat.'+s[3]
                    typemap[key] = dtype
                    shapemap[key] = dshape
                elif s[0] == '#':
                    break

            s = f.read().split('@')
            for key in propmap.keys():
                data = s[key].split('\n')[1:]
                data = ' '.join(data)
                arr = sp.fromstring(data, dtype=typemap[key], sep=' ')
                arr = sp.reshape(arr, newshape=shapemap[key])
                net[propmap[key]] = arr
            # End file parsing

        net['pore.coords'] = net['pore.VertexCoordinates']
        net['throat.conns'] = sp.sort(net['throat.EdgeConnectivity'], axis=1)

        if network is None:
            network = GenericNetwork()
        network = cls._update_network(network=network, net=net)

        return network.project
