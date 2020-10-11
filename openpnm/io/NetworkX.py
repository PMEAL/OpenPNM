import numpy as np
from openpnm.io import GenericIO
from openpnm.utils import logging
from openpnm.network import GenericNetwork
logger = logging.getLogger(__name__)


class NetworkX(GenericIO):
    r"""
    NetworkX is a common tool for dealing with network structures

    This class actually handles dealing with live NetworkX objects rather than
    saved files. NetworkX has a number of export formats, so converting
    OpenPNM data to a NetworkX object, then using NetworkX to exprt might be
    an option for file formats not supported by OpenPNM yet.

    Notes
    -----
    1. Each node in a NetworkX object (i.e. ``G``) can be assigned properties
    using syntax like ``G.node[n]['diameter'] = 0.5`` where ``n`` is the
    node number.  There is no need to precede the property name with any
    indication that it is pore data such as \'pore\_\'.  OpenPNM will prepend
    \'pore.\' to each property name.

    2. Since \'pore.coords\' is so central to OpenPNM it should be specified
    in the NetworkX object as \'coords\', and the [X, Y, Z] coordinates of
    each node should be a 1x3 list.

    3. Edges in a NetworkX object are accessed using the index numbers of the
    two nodes it connects, such as ``G.adj[2][3]['length'] = 0.1``
    indicating the edge that connects nodes 2 and 3.  There is no need to
    precede the property name with any indication that it is throat data such
    as \'throat\_\'.  OpenPNM will prepend \'throat.\' to each property name.

    4. The \'throat.conns\' property is essential to OpenPNM, but this does NOT
    need to be specified explicitly as a property in NetworkX.  The
    connectivity is embedded into the network representation and is extracted
    by OpenPNM.

    """
    @classmethod
    def from_networkx(cls, *args, **kwargs):
        r"""
        This method is being deprecated.  Use ``import_data`` instead.
        """
        return cls.import_data(*args, **kwargs)

    @classmethod
    def import_data(cls, G, project=None):
        r"""
        Add data to an OpenPNM Network from a undirected NetworkX graph object.

        Parameters
        ----------
        G : networkx.classes.graph.Graph Object
            The NetworkX graph. G should be undirected. The numbering of nodes
            should be numeric (int's), zero-based and should not contain any
            gaps, i.e. ``G.nodes() = [0,1,3,4,5]`` is not allowed and should be
            mapped to ``G.nodes() = [0,1,2,3,4]``.

        project : OpenPNM Project object
            A GenericNetwork is created and added to the specified Project.
            If no Project is supplied then one will be created and returned.

        Returns
        -------
        An OpenPNM Project containing a GenericNetwork with all the data from
        the NetworkX object.

        """
        import networkx as nx

        net = {}

        # Ensure G is an undirected networkX graph with numerically numbered
        # nodes for which numbering starts at 0 and does not contain any gaps
        if not isinstance(G, nx.Graph):
            raise Exception('Provided object is not a NetworkX graph.')
        if nx.is_directed(G):
            raise Exception('Provided graph is directed. Convert to undirected graph.')
        if not all(isinstance(n, int) for n in G.nodes()):
            raise Exception('Node numbering is not numeric. Convert to int.')
        if min(G.nodes()) != 0:
            raise Exception('Node numbering does not start at zero.')
        if max(G.nodes()) + 1 != len(G.nodes()):
            raise Exception('Node numbering contains gaps. Map nodes to remove gaps.')

        # Parsing node data
        Np = len(G)
        net.update({'pore.all': np.ones((Np,), dtype=bool)})
        for n, props in G.nodes(data=True):
            for item in props.keys():
                val = props[item]
                dtype = type(val)
                # Remove prepended pore. and pore_ if present
                for b in ['pore.', 'pore_']:
                    item = item.replace(b, '')
                # Create arrays for subsequent indexing, if not present already
                if 'pore.'+item not in net.keys():
                    if dtype == str:  # handle strings of arbitrary length
                        net['pore.'+item] = np.ndarray((Np,), dtype='object')
                    elif dtype is list:
                        dtype = type(val[0])
                        if dtype == str:
                            dtype = 'object'
                        cols = len(val)
                        net['pore.'+item] = np.ndarray((Np, cols), dtype=dtype)
                    else:
                        net['pore.'+item] = np.ndarray((Np,), dtype=dtype)
                net['pore.'+item][n] = val

        # Parsing edge data
        # Deal with conns explicitly
        try:
            conns = list(G.edges)   # NetworkX V2
        except Exception:
            conns = G.edges()       # NetworkX V1
        conns.sort()

        # Add conns to Network
        Nt = len(conns)
        net.update({'throat.all': np.ones(Nt, dtype=bool)})
        net.update({'throat.conns': np.array(conns)})

        # Scan through each edge and extract all its properties
        i = 0
        for t in conns:
            props = G[t[0]][t[1]]
            for item in props:
                val = props[item]
                dtype = type(val)
                # Remove prepended throat. and throat_ if present
                for b in ['throat.', 'throat_']:
                    item = item.replace(b, '')
                # Create arrays for subsequent indexing, if not present already
                if 'throat.'+item not in net.keys():
                    if dtype == str:
                        net['throat.'+item] = np.ndarray((Nt,), dtype='object')
                    if dtype is list:
                        dtype = type(val[0])
                        if dtype == str:
                            dtype = 'object'
                        cols = len(val)
                        net['throat.'+item] = np.ndarray((Nt, cols),
                                                         dtype=dtype)
                    else:
                        net['throat.'+item] = np.ndarray((Nt,), dtype=dtype)
                net['throat.'+item][i] = val
            i += 1

        network = GenericNetwork(project=project)
        network = cls._update_network(network=network, net=net)
        return network.project

    @classmethod
    def to_networkx(cls, *args, **kwargs):
        r"""
        This method is being deprecated.  Use ``export_data`` instead.
        """
        return cls.export_data(*args, **kwargs)

    @classmethod
    def export_data(cls, network):
        r"""
        Write OpenPNM Network to a NetworkX object.

        Parameters
        ----------
        network : OpenPNM Network Object
            The OpenPNM Network to be converted to a NetworkX object

        Returns
        -------
        A NetworkX object with all pore/throat properties attached to it
        """
        import networkx as nx

        # Ensure network is an OpenPNM Network object.
        if not isinstance(network, GenericNetwork):
            raise('Provided network is not an OpenPNM Network.')

        G = nx.Graph()

        # Extracting node list and connectivity matrix from Network
        nodes = map(int, network.Ps)
        conns = network['throat.conns']

        # Explicitly add nodes and connectivity matrix
        G.add_nodes_from(nodes)
        G.add_edges_from(conns)

        # Attach Network properties to G
        for prop in network.props(deep=True) + network.labels():
            if 'pore.' in prop:
                if len(network[prop].shape) > 1:
                    val = {i: list(network[prop][i]) for i in network.Ps}
                else:
                    val = {i: network[prop][i] for i in network.Ps}
                nx.set_node_attributes(G, name=prop[5:], values=val)
            if 'throat.' in prop:
                val = {tuple(conn): network[prop][i] for i, conn
                       in enumerate(conns)}
                nx.set_edge_attributes(G, name=prop[7:], values=val)

        return G
