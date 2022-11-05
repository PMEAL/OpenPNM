import logging
import numpy as np
from openpnm.io import _parse_filename
from openpnm.network import Network


logger = logging.getLogger(__name__)


def network_to_pergeos(network, filename=''):
    r"""
    """
    # avoid printing truncated array
    np.set_printoptions(threshold=np.inf)

    # Ensure network has PerGeos' expected properties
    if 'pore.EqRadius' not in network.props():
        try:
            network['pore.EqRadius'] = network['pore.diameter']/2
        except KeyError:
            network['pore.EqRadius'] = np.ones([network.Np, ])

    s = ["# Avizo 3D ASCII 3.0\n\n"]
    s.append("define VERTEX " + str(network.Np) + '\n')
    s.append("define EDGE " + str(network.Nt) + '\n')
    s.append("define POINT " + str(2*network.Nt) + '\n\n')
    s.append("Parameters {\n\tContentType \"HxPoreNetworkModel\"\n}\n\n")

    types = {'b': 'int', 'i': 'int', 'f': 'float'}
    typemap = {}
    namemap = {}
    shapemap = {}
    propmap = {}
    i = 1

    NumEdgePoints = 1
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
        temp = element[1] + " { " + typemap[item] + shapemap[item] + " " \
            + namemap[item] + " } @" + str(i) + '\n'

        if temp.find('EdgeConnectivity') == -1:
            # replaces openpnm tags with the mandatory am file's tags
            if "Conns" in temp:
                temp = temp.replace("Conns", "EdgeConnectivity")
            elif "Coords" in temp:
                temp = temp.replace("Coords", "VertexCoordinates")
            s.append(temp)
        propmap[item] = str(i)
        if "NumEdgePoints" in temp:
            NumEdgePoints = 0
        i += 1

    if NumEdgePoints:
        temp = "EDGE { int NumEdgePoints" + " } @" + str(i) + '\n'
        s.append(temp)
        tempat = "@" + str(i) + '\n'
        i += 1

    # Add POINT data
    s.append("POINT { float[3] EdgePointCoordinates } @" + str(i))
    s.append("\n\n# Data section follows")
    for item in network.keys():
        data = network[item]
        if item != 'throat.EdgeConnectivity':
            s.append('\n\n@' + propmap[item] + '\n')
            if shapemap[item] == '':
                data = np.atleast_2d(data).T
            if typemap[item] == 'float':
                formatter = {'float_kind': lambda x: "%.15E" % x}
            else:
                formatter = None
            if data.dtype == 'bool':
                data = data.astype(int)
            d = np.array2string(data, formatter=formatter)
            s.append(d.replace('[', '').replace(']', '').replace('\n ', '\n'))

    # Add POINT data
    s.append('\n\n@' + str(i) + '\n')
    formatter = {'float_kind': lambda x: "%.15E" % x}

    conns = network['throat.conns']
    d = np.array2string(network['pore.coords'][conns], formatter=formatter)
    for r in (('[', ''), (']', ''), ('\n\n', '\n'), ('\n  ', '\n'),
              ('\n ', '\n')):
        d = d.replace(*r)
    d += '\n'
    s.append(d)

    # Add NumEdgePoints
    if NumEdgePoints:
        s.append('\n\n' + tempat)
        s.append(''.join(['2' + '\n']*network.Nt))

    # Write to file
    if filename == '':
        filename = network.name
    fname = _parse_filename(filename=filename, ext='am')
    with open(fname, 'w') as f:
        f.write(''.join(s))


def network_from_pergeos(filename):
    r"""

    Notes
    -----
    PerGeos is the format used by the Avizo software. See `here for more
    details <https://cases.pergeos.com/>`_.
    """
    net = {}

    # ---------------------------------------------------------------------
    # Parse the link1 file
    filename = _parse_filename(filename=filename, ext='am')
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
                temp = np.zeros(dshape, dtype=dtype)
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
                temp = np.zeros(dshape, dtype=dtype)
                net['throat.'+s[3]] = temp
                key = int(s[-1].replace('@', ''))
                propmap[key] = 'throat.'+s[3]
                typemap[key] = dtype
                shapemap[key] = dshape
            elif s[0] == '#':
                break

        s = f.read().split('@')
        for key in propmap.keys():
            if key in s:
                data = s[key].split('\n')[1:]
                data = ' '.join(data)
                arr = np.fromstring(data, dtype=typemap[key], sep=' ')
                arr = np.reshape(arr, newshape=shapemap[key])
                net[propmap[key]] = arr
        # End file parsing

    net['pore.coords'] = net['pore.VertexCoordinates']
    net['throat.conns'] = np.sort(net['throat.EdgeConnectivity'], axis=1)

    network = Network()
    network.update(net)

    return network
