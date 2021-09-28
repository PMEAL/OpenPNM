import numpy as np
import scipy as sp
import openpnm as op
from networkx import complete_graph, random_layout
from networkx import set_node_attributes, set_edge_attributes


class NetworkXTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.settings['local_data'] = True
        self.net = op.network.Cubic(shape=[2, 2, 2])
        Ps = [0, 1, 2, 3]
        Ts = self.net.find_neighbor_throats(pores=Ps)
        self.geo_1 = op.geometry.GenericGeometry(network=self.net,
                                                 pores=Ps, throats=Ts)
        self.geo_1['pore.boo'] = 1
        self.geo_1['throat.boo'] = 1
        Ps = [4, 5, 6, 7]
        Ts = self.net.find_neighbor_throats(pores=Ps, mode='xnor')
        self.geo_2 = op.geometry.GenericGeometry(network=self.net,
                                                 pores=Ps, throats=Ts)
        self.geo_2['pore.boo'] = 1
        self.geo_2['throat.boo'] = 1

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_from_networkx(self):
        G = complete_graph(10)
        pos = random_layout(G, dim=3)
        val = {n: list(pos[n]) for n in pos}
        set_node_attributes(G, name='coords', values=val)
        set_node_attributes(G, name='area', values=1.123)
        set_node_attributes(G, name='diameter', values=1.123)
        set_edge_attributes(G, name='length', values=1.123)
        set_edge_attributes(G, name='perimeter', values=1.123)
        project = op.io.NetworkX.from_networkx(G=G)
        assert len(project) == 1
        num_nodes = len(G.nodes())
        num_edges = len(G.edges())
        net = project.network
        assert net.Np == num_nodes
        assert net.Nt == num_edges
        assert np.shape(net['pore.coords']) == (num_nodes, 3)
        assert np.shape(net['throat.conns']) == (num_edges, 2)
        a = {'pore.area', 'pore.diameter', 'throat.length', 'throat.perimeter'}
        assert a.issubset(net.props())

    def test_save_and_load_networkx_no_phases(self):
        G = op.io.NetworkX.export_data(network=self.net)
        project = op.io.NetworkX.import_data(G)
        assert len(project) == 1
        net = project.network
        assert net.Np == 8
        assert net.Nt == 12
        assert np.shape(net['pore.coords']) == (8, 3)
        assert np.shape(net['throat.conns']) == (12, 2)


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = NetworkXTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
