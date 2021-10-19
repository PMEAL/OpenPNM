import numpy as np
import openpnm as op


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

    def test_nearest_and_furthest_neighbor(self):
        pn = op.network.Cubic(shape=[3, 3, 3], spacing=[1, 2, 3])
        pn.add_model(propname='pore.L',
                     model=op.models.topology.distance_to_nearest_neighbor)
        assert np.all(pn['pore.L'][:5] == np.array([1, 1, 1, 1, 1]))
        pn.add_model(propname='pore.L',
                     model=op.models.topology.distance_to_furthest_neighbor)
        assert np.all(pn['pore.L'][:5] == np.array([3, 3, 3, 3, 3]))

    def test_nearest_pore(self):
        pn = op.network.Cubic(shape=[3, 3, 3], spacing=1)
        pn.add_model(propname='pore.nearby',
                     model=op.models.topology.distance_to_nearest_pore)
        assert pn['pore.nearby'][1] == 1.0
        pn['pore.coords'][1] = pn['pore.coords'][0]
        pn.regenerate_models('pore.nearby')
        assert pn['pore.nearby'][0] == 0.0
        assert pn['pore.nearby'][1] == 0.0

    def test_find_coincident_pores(self):
        pn = op.network.Cubic(shape=[3, 3, 3], spacing=1)
        pn.add_model(propname='pore.nearby',
                     model=op.models.topology.find_coincident_pores)
        assert pn['pore.nearby'][0] == []
        pn['pore.coords'][1] = pn['pore.coords'][0]
        pn.regenerate_models('pore.nearby')
        assert pn['pore.nearby'][0] == [1]
        assert pn['pore.nearby'][1] == [0]
        assert pn['pore.nearby'][2] == []

    def test_count_coincident_pores(self):
        pn = op.network.Cubic(shape=[3, 3, 3], spacing=1)
        pn.add_model(propname='pore.nearby',
                     model=op.models.topology.count_coincident_pores)
        assert pn['pore.nearby'][0] == 0
        pn['pore.coords'][1] = pn['pore.coords'][0]
        pn.regenerate_models('pore.nearby')
        assert pn['pore.nearby'][0] == 1
        assert pn['pore.nearby'][1] == 1
        assert pn['pore.nearby'][2] == 0

    def test_cluster_number(self):
        pn = op.network.Cubic(shape=[6, 1, 1])
        op.topotools.trim(network=pn, throats=[0, 3])
        pn.add_model(propname='pore.cluster',
                     model=op.models.topology.cluster_number)
        assert np.all(pn['pore.cluster'] == [0, 1, 1, 1, 2, 2])

    def test_cluster_size(self):
        pn = op.network.Cubic(shape=[6, 1, 1])
        op.topotools.trim(network=pn, throats=[0, 3])
        pn.add_model(propname='pore.cluster_size',
                     model=op.models.topology.cluster_size)
        assert np.all(pn['pore.cluster_size'] == [1, 3, 3, 3, 2, 2])
        pn.add_model(propname='pore.cluster_number',
                     model=op.models.topology.cluster_number)
        pn.regenerate_models()
        assert np.all(pn['pore.cluster_size'] == [1, 3, 3, 3, 2, 2])

    def test_duplicate_throats(self):
        pn = op.network.Cubic(shape=[6, 1, 1])
        op.topotools.extend(network=pn, throat_conns=[[0, 1]])
        pn.add_model(propname='throat.dupes',
                     model=op.models.topology.duplicate_throats)
        assert pn['throat.dupes'].sum() == 1

    def test_headless_throats(self):
        pn = op.network.Cubic(shape=[6, 1, 1])
        pn['throat.conns'][0, :] = [0, 110]
        pn.add_model(propname='throat.headless',
                     model=op.models.topology.headless_throats)
        assert pn['throat.headless'].sum() == 1

    def test_isolated_pores(self):
        pn = op.network.Cubic(shape=[6, 1, 1])
        op.topotools.trim(network=pn, throats=[0])
        pn.add_model(propname='pore.isolated',
                     model=op.models.topology.isolated_pores)
        assert pn['pore.isolated'].sum() == 1

    def test_looped_throats(self):
        pn = op.network.Cubic(shape=[6, 1, 1])
        op.topotools.extend(network=pn, throat_conns=[[1, 1]])
        pn.add_model(propname='throat.looped',
                     model=op.models.topology.looped_throats)
        assert pn['throat.looped'].sum() == 1

    def test_reversed_throats(self):
        pn = op.network.Cubic(shape=[6, 1, 1])
        pn['throat.conns'][0, :] = [1, 0]
        pn.add_model(propname='throat.reversed',
                     model=op.models.topology.reversed_throats)
        assert pn['throat.reversed'].sum() == 1


if __name__ == '__main__':

    t = TopologyModelsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
