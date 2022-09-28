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

    def test_ispercolating_bond(self):
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

    def test_remove_isolated_clusters(self):
        pass

    def test_site_percolation(self):
        pass

    def test_bond_percolation(self):
        pass

    def test_trim_disconnected_clusters(self):
        pass

    def test_find_clusters(self):
        pass



if __name__ == '__main__':

    t = PerctoolsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: ' + item)
            t.__getattribute__(item)()
