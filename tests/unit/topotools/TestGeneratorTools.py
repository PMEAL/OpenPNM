import pytest
import numpy as np
import openpnm as op
import matplotlib.pyplot as plt
from openpnm import topotools


class GeneratorToolTest:

    def setup_class(self):
        self.ws = op.Workspace()

    def test_trim_verts(self):
        net = op.topotools.generators.cubic(shape=[3, 3, 3])
        net['vert.test'] = np.ones(net['vert.coords'].shape[0])
        net['edge.test'] = np.ones(net['edge.conns'].shape[0])
        assert net['vert.test'].sum() == 27
        assert net['edge.test'].sum() == 54
        net = op.topotools.generators.tools.trim(network=net,
                                                 vert_ids=[4])
        assert net['vert.coords'].shape[0] == 26
        assert net['edge.conns'].max() == 25
        assert net['vert.test'].sum() == 26
        assert net['edge.test'].sum() == 49

    def test_trim_edges(self):
        net = op.topotools.generators.cubic(shape=[3, 3, 3])
        net['vert.test'] = np.ones(net['vert.coords'].shape[0])
        net['edge.test'] = np.ones(net['edge.conns'].shape[0])
        assert net['vert.test'].sum() == 27
        assert net['edge.test'].sum() == 54
        net = op.topotools.generators.tools.trim(network=net,
                                                 edge_ids=[4])
        assert net['vert.coords'].shape[0] == 27
        assert net['edge.conns'].shape[0] == 53
        assert net['edge.conns'].max() == 26
        assert net['vert.test'].sum() == 27
        assert net['edge.test'].sum() == 53

    def test_trim_edges_and_verts(self):
        net = op.topotools.generators.cubic(shape=[3, 3, 3])
        with pytest.raises(Exception):
            net = op.topotools.generators.tools.trim(network=net,
                                                     edge_ids=[4],
                                                     vert_ids=[4])

    def test_join(self):
        net1 = op.topotools.generators.cubic(shape=[3, 3, 3])
        net2 = op.topotools.generators.cubic(shape=[3, 3, 3])
        net2['vert.coords'] += 3
        net1['edge.test'] = np.ones(net1['edge.conns'].shape[0])
        net2['edge.test'] = np.ones(net2['edge.conns'].shape[0])
        net1['vert.test'] = np.ones((net1['vert.coords'].shape[0], 2))
        net2['vert.test'] = np.ones((net2['vert.coords'].shape[0], 2))
        net3 = op.topotools.generators.tools.join(net1, net2)
        assert net3['edge.conns'].max() == 53
        assert net3['vert.coords'].shape[0] == 54
        assert net3['edge.conns'].shape[0] == 108
        assert np.all(net3['vert.test'].shape == np.array((54, 2)))
        assert np.all(net3['edge.test'].shape == np.array([108,]))

    def test_get_spacing(self):
        net = op.topotools.generators.cubic(shape=[3, 3, 3], spacing=[1, 2, 3])
        spc = op.topotools.generators.tools.get_spacing(net)
        assert np.all(np.array([1, 2, 3]) == spc)

    def test_get_shape(self):
        net = op.topotools.generators.cubic(shape=[3, 4, 5])
        shp = op.topotools.generators.tools.get_shape(net)
        assert np.all(np.array([3, 4, 5]) == shp)



if __name__ == '__main__':

    t = GeneratorToolTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
