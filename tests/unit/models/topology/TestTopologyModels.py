import numpy as np
import openpnm as op
from numpy.testing import assert_allclose


class TopologyModelsTest:

    def setup_class(self):
        pass

    def test_coordination_number(self):
        pn = op.network.Cubic(shape=[3, 3, 3])
        pn.add_model(propname='pore.z',
                     model=op.models.topology.coordination_number)
        assert np.all(pn['pore.z'][:5] == np.array([3, 4, 3, 4, 5]))

    def test_pore_to_pore_distance(self):
        pn = op.network.Cubic(shape=[3, 3, 3], spacing=[1, 2, 3])
        pn.add_model(propname='throat.L',
                     model=op.models.topology.pore_to_pore_distance)
        assert np.all(pn['throat.L'][[0, 18, 50]] == np.array([3, 2, 1]))

    def test_nearest_neighbor(self):
        pn = op.network.Cubic(shape=[3, 3, 3], spacing=[1, 2, 3])
        pn.add_model(propname='pore.L',
                     model=op.models.topology.nearest_neighbor_distance)
        assert np.all(pn['pore.L'][:5] == np.array([1, 1, 1, 1, 1]))
        pn.add_model(propname='pore.L',
                     model=op.models.topology.furthest_neighbor_distance)
        assert np.all(pn['pore.L'][:5] == np.array([3, 3, 3, 3, 3]))


if __name__ == '__main__':

    t = TopologyModelsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
