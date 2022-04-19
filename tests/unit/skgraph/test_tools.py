import pytest
import numpy as np
from numpy.testing import assert_allclose
from openpnm._skgraph import tools
from openpnm._skgraph.generators import cubic


class SKGRToolsTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_dict_to_am_undirected(self):
        g = cubic(shape=[3, 2, 1])
        am = tools.dict_to_am(g)
        # Make sure am is symmetrical but edges have the same order
        assert np.all(am.row == np.hstack((g['edge.conns'][:, 0], g['edge.conns'][:, 1])))
        assert np.all(am.col == np.hstack((g['edge.conns'][:, 1], g['edge.conns'][:, 0])))
        assert_allclose(np.linalg.norm(am.todense()), 3.7416573867739413)

    def test_dict_to_am_directed(self):
        g = cubic(shape=[3, 2, 1])
        g['edge.conns'][1, :] = [1, 0]
        am = tools.dict_to_am(g)
        # Make sure edges in am are untouched
        assert np.all(am.row == g['edge.conns'][:, 0])
        assert np.all(am.col == g['edge.conns'][:, 1])
        assert_allclose(np.linalg.norm(am.todense()), 2.6457513110645907)

    def test_dict_to_am_undirected_w_weights(self):
        g = cubic(shape=[3, 2, 1])
        Ts = np.arange(g['edge.conns'].shape[0])
        am = tools.dict_to_am(g, weights=Ts)
        assert np.all(am.data == np.hstack((Ts, Ts)))

    def test_dict_to_am_directed_w_weights(self):
        g = cubic(shape=[3, 2, 1])
        g['edge.conns'][1, :] = [1, 0]
        Ts = np.arange(g['edge.conns'].shape[0])
        am = tools.dict_to_am(g, weights=Ts)
        assert np.all(am.data == Ts)

    def test_dict_to_am_w_dupes(self):
        g = cubic(shape=[3, 2, 1])
        g['edge.conns'][1, :] = [0, 1]
        # This does NOT raise an exception since checking for this would add
        # too much overhead
        am = tools.dict_to_am(g)
        assert len(am.data == 14)
        assert am.data.max() == 1
        # If we remove duplicates now, am is changed
        am.sum_duplicates()
        assert len(am.data == 12)
        assert am.data.max() == 2

    def test_dict_to_im_undirected(self):
        g = cubic(shape=[3, 2, 1])
        im = tools.dict_to_im(g)
        for i in range(im.shape[0]):
            a = np.where(im.getcol(i).todense())[0]
            b = g['edge.conns'][i, :]
            assert np.all(a == b)

    def test_dict_to_im_directed(self):
        g = cubic(shape=[3, 2, 1])
        g['edge.conns'] = [3, 2]
        with pytest.raises(Exception):
            _ = tools.dict_to_im(g)

    def test_cart2cyl_and_back(self):
        x, y, z = np.random.rand(10, 3).T
        r, q, z1 = tools.cart2cyl(x, y, z)
        x2, y2, z2 = tools.cyl2cart(r, q, z1)
        assert_allclose(x, x2)
        assert_allclose(y, y2)
        assert_allclose(z, z2)

    def test_cart2cyl_polar(self):
        x, y = np.random.rand(10, 2).T
        r, q, z = tools.cart2cyl(x, y)
        assert z.sum() == 0
        x, y, z = tools.cyl2cart(r, q)
        assert z.sum() == 0

    def test_cart2sph_and_back(self):
        x, y, z = np.random.rand(10, 3).T
        r, q, phi = tools.cart2sph(x, y, z)
        x2, y2, z2 = tools.sph2cart(r, q, phi)
        assert_allclose(x, x2)
        assert_allclose(y, y2)
        assert_allclose(z, z2)

    def test_find_coincident_nodes(self):
        g = cubic(shape=[10, 10, 10])
        hits = tools.find_coincident_nodes(g)
        assert hits == []
        g['node.coords'][1, :] = g['node.coords'][0, :]
        g['node.coords'][2, :] = g['node.coords'][0, :]
        g['node.coords'][4, :] = g['node.coords'][0, :]
        hits = tools.find_coincident_nodes(g)
        assert len(hits) == 1
        assert np.all(hits[0] == [0, 1, 2, 4])
        g['node.coords'][10, :] = g['node.coords'][20, :]
        hits = tools.find_coincident_nodes(g)
        assert len(hits) == 2
        assert np.all(hits[1] == [10, 20])

    def test_dimensionality(self):
        g = cubic(shape=[3, 1, 1])
        assert np.all(tools.dimensionality(g) == [True, False, False])
        g = cubic(shape=[1, 3, 1])
        assert np.all(tools.dimensionality(g) == [False, True, False])
        g = cubic(shape=[1, 1, 3])
        assert np.all(tools.dimensionality(g) == [False, False, True])
        g = cubic(shape=[3, 3, 1])
        assert np.all(tools.dimensionality(g) == [True, True, False])
        g = cubic(shape=[3, 1, 3])
        assert np.all(tools.dimensionality(g) == [True, False, True])
        g = cubic(shape=[1, 3, 3])
        assert np.all(tools.dimensionality(g) == [False, True, True])
        g = cubic(shape=[3, 3, 3])
        assert np.all(tools.dimensionality(g) == [True, True, True])

    def test_get_prefixes(self):
        g = cubic(shape=[3, 1, 1], node_prefix='foo', edge_prefix='bar')
        node_prefix = tools.get_node_prefix(g)
        edge_prefix = tools.get_edge_prefix(g)
        assert node_prefix == 'foo'
        assert edge_prefix == 'bar'
        # The following should be made to throw an exception, but is a tiny
        # edge case that would add a lot of overhead since it would mean
        # scanning the entire keys list rather than just stopping at the first
        # occurance of 'coords' (or 'conns').
        # g['test.coords'] = 1
        # node_prefix = tools.get_node_prefix(g)
        g['foo.coords2'] = 2
        node_prefix = tools.get_node_prefix(g)
        assert node_prefix == 'foo'
        g = tools.change_prefix(g, 'bar', 'baz')
        edge_prefix = tools.get_edge_prefix(g)
        assert edge_prefix == 'baz'

    def test_get_spacing(self):
        d = cubic(shape=[3, 3, 3])
        s = tools.get_cubic_spacing(d)
        assert np.all(s == [1, 1, 1])
        d = cubic(shape=[3, 3, 3], spacing=[1, 2, 3])
        s = tools.get_cubic_spacing(d)
        assert np.all(s == [1, 2, 3])
        d = cubic(shape=[5, 4, 3], spacing=[1, 2, 3])
        s = tools.get_cubic_spacing(d)
        assert np.all(s == [1, 2, 3])

    def test_get_shape(self):
        d = cubic(shape=[3, 3, 3])
        s = tools.get_cubic_shape(d)
        assert np.all(s == [3, 3, 3])
        d = cubic(shape=[3, 3, 3], spacing=[1, 2, 3])
        s = tools.get_cubic_shape(d)
        assert np.all(s == [3, 3, 3])
        d = cubic(shape=[5, 4, 3], spacing=[1, 2, 3])
        s = tools.get_cubic_shape(d)
        assert np.all(s == [5, 4, 3])


if __name__ == '__main__':
    t = SKGRToolsTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
