import pytest
import numpy as np
import openpnm as op
from numpy.testing import assert_allclose
from openpnm import topotools


class PerctoolsTest:

    def setup_class(self):
        self.ws = op.Workspace()

    def teardown_class(self):
        self.ws.clear()

    def test_find_path(self):
        pn = op.network.Demo(shape=[4, 4, 1])
        nodes, edges = op.topotools.find_path(
            network=pn,
            pore_pairs=[[0, 15]]).values()
        assert len(nodes) == 1
        assert len(edges) == 1
        assert np.all(nodes[0] == [0, 4, 8, 9, 13, 14, 15])
        assert np.all(edges[0] == [12, 16, 6, 21, 10, 11])
        nodes, edges = op.topotools.find_path(
            network=pn,
            pore_pairs=[[0, 15], [1, 15]]).values()
        assert len(nodes) == 2
        assert len(edges) == 2
        assert np.all(nodes[0] == [0, 4, 8, 9, 13, 14, 15])
        assert np.all(edges[0] == [12, 16, 6, 21, 10, 11])

    def test_ispercolating(self):
        pn = op.network.Demo(shape=[4, 4, 1])
        flag = op.topotools.ispercolating(
            network=pn,
            inlets=pn.pores('left'),
            outlets=pn.pores('right'))
        assert flag
        op.topotools.trim(network=pn, pores=[1, 5, 9, 13])
        flag = op.topotools.ispercolating(
            network=pn,
            inlets=pn.pores('left'),
            outlets=pn.pores('right'))
        assert flag
        op.topotools.trim(network=pn, pores=[3, 4, 5])
        flag = op.topotools.ispercolating(
            network=pn,
            inlets=pn.pores('left'),
            outlets=pn.pores('right'))
        assert not flag

    def test_find_isolated_clusters_pore_mask(self):
        pn = op.network.Demo(shape=[4, 4, 1])
        plabels = np.zeros_like(pn.Ps, dtype=bool)
        np.random.seed(0)
        hits = np.random.randint(0, pn.Np, 8)
        plabels[hits] = True
        labels = op.topotools.find_isolated_clusters(
            network=pn,
            mask=plabels,
            inlets=pn.pores('left')
        )
        assert np.all(labels == [5, 12])
        # ax = op.visualization.plot_coordinates(pn, pores=labels,
        #                                        s=50, c='r')
        # ax = op.visualization.plot_coordinates(pn, pores=plabels,
        #                                        s=50, marker='x', c='b', ax=ax)
        # ax = op.visualization.plot_coordinates(pn, pores=pn.pores('left'),
        #                                        s=50, marker='x', c='g', ax=ax)

    def test_find_isolated_clusters_throat_mask(self):
        pn = op.network.Demo(shape=[4, 4, 1])
        mask = np.zeros_like(pn.Ts, dtype=bool)
        np.random.seed(0)
        hits = np.random.randint(0, pn.Nt, 8)
        mask[hits] = True
        inlets = pn.pores('left')
        labels = op.topotools.find_isolated_clusters(
            network=pn,
            mask=mask,
            inlets=inlets
        )
        # ax = op.visualization.plot_connections(pn, throats=mask)
        # ax = op.visualization.plot_coordinates(pn, pores=labels,
        #                                         s=50, c='r', ax=ax)
        # ax = op.visualization.plot_coordinates(pn, pores=plabels,
        #                                         s=50, marker='x', c='b', ax=ax)
        # ax = op.visualization.plot_coordinates(pn, pores=pn.pores('left'),
        #                                         s=50, marker='x', c='g', ax=ax)

    def test_site_percolation(self):
        pass

    def test_bond_percolation(self):
        pass

    def test_trim_disconnected_clusters(self):
        pass

    def test_find_clusters_sites(self):
        net = op.network.Cubic(shape=[10, 10, 1])
        net['pore.seed'] = np.random.rand(net.Np)
        net['throat.seed'] = np.random.rand(net.Nt)
        clusters = topotools.find_clusters(network=net,
                                           mask=net['pore.seed'] < 0.5)
        assert len(clusters.pore_labels) == net.Np
        assert len(clusters.throat_labels) == net.Nt
        # Ensure neighboring pores have same label, unless one is -1
        L = clusters.pore_labels[net.conns]
        hits = np.all(L >= 0, axis=1)
        assert np.all(L[:, 0][hits] == L[:, 1][hits])

    def test_find_clusters_bonds(self):
        net = op.network.Cubic(shape=[10, 10, 1])
        net['pore.seed'] = np.random.rand(net.Np)
        net['throat.seed'] = np.random.rand(net.Nt)
        clusters = topotools.find_clusters(network=net,
                                           mask=net['throat.seed'] < 0.5)
        assert len(clusters.pore_labels) == net.Np
        assert len(clusters.throat_labels) == net.Nt
        # Ensure neighboring pores have same label, if throat is invaded
        L = clusters.pore_labels[net.conns]
        hits = net['throat.seed'] < 0.5
        assert np.all(L[:, 0][hits] == L[:, 1][hits])


if __name__ == '__main__':

    t = PerctoolsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
