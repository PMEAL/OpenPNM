import numpy as np
import scipy as sp
import openpnm as op


class GenericNetworkTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[10, 10, 10])

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_find_connected_pores_numeric_not_flattend(self):
        a = self.net.find_connected_pores(throats=[0, 1])
        assert np.all(a.flatten() == [0, 1, 1, 2])

    def test_find_connected_pores_numeric_flattend(self):
        a = self.net.find_connected_pores(throats=[0, 1], flatten=True)
        assert np.all(a == [0, 1, 2])

    def test_find_connected_pores_boolean_flattend(self):
        Tind = np.zeros((self.net.Nt,), dtype=bool)
        Tind[[0, 1]] = True
        a = self.net.find_connected_pores(throats=Tind, flatten=True)
        assert np.all(a == [0, 1, 2])

    def test_find_connected_pores_empty_flattend(self):
        a = self.net.find_connected_pores(throats=[], flatten=True)
        assert np.shape(a) == (0, )

    def test_find_neighbor_pores_numeric(self):
        a = self.net.find_neighbor_pores(pores=[])
        assert np.size(a) == 0

    def test_find_neighbor_pores_boolean(self):
        Pind = np.zeros((self.net.Np,), dtype=bool)
        Pind[[0, 1]] = True
        a = self.net.find_neighbor_pores(pores=Pind)
        assert np.all(a == [2, 10, 11, 100, 101])

    def test_find_neighbor_pores_numeric_union(self):
        a = self.net.find_neighbor_pores(pores=[0, 2], mode='union')
        assert np.all(a == [1, 3, 10, 12, 100, 102])

    def test_find_neighbor_pores_numeric_intersection(self):
        a = self.net.find_neighbor_pores(pores=[0, 2], mode='xnor')
        assert np.all(a == [1])

    def test_find_neighbor_pores_numeric_exclusive_or(self):
        a = self.net.find_neighbor_pores(pores=[0, 2], mode='exclusive_or')
        assert np.all(a == [3, 10, 12, 100, 102])

    def test_find_neighbor_pores_numeric_union_include_input(self):
        a = self.net.find_neighbor_pores(pores=[0, 2], mode='or',
                                         include_input=True)
        assert np.all(a == [1, 3, 10, 12, 100, 102])
        a = self.net.find_neighbor_pores(pores=[0, 1], mode='or',
                                         include_input=True)
        assert np.all(a == [0, 1, 2, 10, 11, 100, 101])

    def test_find_neighbor_pores_numeric_intersection_include_input(self):
        a = self.net.find_neighbor_pores(pores=[0, 2], mode='and',
                                         include_input=True)
        assert np.all(a == [1])
        a = self.net.find_neighbor_pores(pores=[0, 1], mode='and',
                                         include_input=True)
        assert np.all(a == [])

    def test_find_neighbor_pores_numeric_intersection_exclude_input(self):
        a = self.net.find_neighbor_pores(pores=[0, 2], mode='and',
                                         include_input=False)
        assert np.all(a == [1])

    def test_find_neighbor_pores_numeric_exclusive_or_include_input(self):
        a = self.net.find_neighbor_pores(pores=[0, 2], mode='exclusive_or',
                                         include_input=True)
        assert np.all(a == [3, 10, 12, 100, 102])
        a = self.net.find_neighbor_pores(pores=[0, 1], mode='exclusive_or',
                                         include_input=True)
        assert np.all(a == [0, 1, 2, 10, 11, 100, 101])

    def test_find_neighbor_throats_on_pores_wo_throats(self):
        net = op.network.Cubic(shape=[10, 10, 1])
        ts = net.find_neighbor_throats(pores=net.Ps[-1])
        op.topotools.trim(net, throats=ts)
        ts2 = net.find_neighbor_throats(pores=99)
        assert ts2.size == 0

    def test_find_neighbor_throats_empty(self):
        a = self.net.find_neighbor_throats(pores=[])
        assert np.size(a) == 0

    def test_find_neighbor_throats_boolean(self):
        Pind = np.zeros((self.net.Np,), dtype=bool)
        Pind[[0, 1]] = True
        a = self.net.find_neighbor_throats(pores=Pind)
        assert np.all(a == [0, 1, 900, 901, 1800, 1801])

    def test_find_neighbor_throats_numeric_union(self):
        a = self.net.find_neighbor_throats(pores=[0, 2], mode='union')
        assert np.all(a == [0, 1, 2, 900, 902, 1800, 1802])

    def test_find_neighbor_throats_numeric_intersection(self):
        a = self.net.find_neighbor_throats(pores=[0, 2], mode='xnor')
        assert np.size(a) == 0

    def test_find_neighbor_throats_numeric_exclusive_or(self):
        a = self.net.find_neighbor_throats(pores=[0, 2],
                                           mode='exclusive_or')
        assert np.all(a == [0, 1, 2, 900, 902, 1800, 1802])

    def test_num_neighbors_empty(self):
        a = self.net.num_neighbors(pores=[])
        assert np.size(a) == 0

    def test_num_neighbors_pores_flattened(self):
        a = self.net.num_neighbors(pores=0, flatten=True)
        assert a == 3
        assert isinstance(a, int)
        a = self.net.num_neighbors(pores=[0, 2], flatten=True)
        assert a == 6
        assert isinstance(a, int)

    def test_num_neighbors_pores_with_modes(self):
        a = self.net.num_neighbors(pores=[0, 2], mode='union', flatten=True)
        assert a == 6
        a = self.net.num_neighbors(pores=[0, 2], mode='xnor',
                                   flatten=True)
        assert a == 1
        a = self.net.num_neighbors(pores=[0, 2], mode='exclusive_or',
                                   flatten=True)
        assert a == 5

    def test_num_neighbors_pores_not_flattened(self):
        a = self.net.num_neighbors(pores=[0, 2], flatten=False)
        assert np.all(a == [3, 4])
        a = self.net.num_neighbors(pores=0, flatten=False)
        assert np.all(a == [3])
        assert isinstance(a, np.ndarray)

    def test_find_nearby_pores_distance_1(self):
        a = self.net.find_nearby_pores(pores=[0, 1], r=1, flatten=False,
                                       include_input=True)
        b = self.net.find_neighbor_pores(pores=[0, 1], flatten=False,
                                         include_input=True)
        assert np.all([np.all(a[i] == b[i]) for i in range(0, len(a))])

    def test_find_nearby_pores_distance_2(self):
        a = self.net.find_nearby_pores(pores=[0, 1], r=2)
        assert np.all([np.size(a[i]) for i in [0, 1]] == [9, 13])

    def test_find_nearby_pores_distance_0(self):
        a = self.net.find_nearby_pores(pores=[0, 1], r=1e-9, flatten=False)
        assert np.shape(a) == (2, 0)
        a = self.net.find_nearby_pores(pores=[0, 1], r=1e-9, flatten=True)
        assert a.shape == (0,)

    def test_find_nearby_pores_distance_1_flattened(self):
        a = self.net.find_nearby_pores(pores=[0, 1], r=1, flatten=True)
        b = self.net.find_neighbor_pores(pores=[0, 1])
        assert np.all(a == b)

    def test_find_nearby_pores_distance_2_flattened(self):
        a = self.net.find_nearby_pores(pores=[0, 1], r=2, flatten=True)
        assert np.size(a) == 15

    def test_find_nearby_pores_distance_2_flattened_include_input(self):
        a = self.net.find_nearby_pores(pores=[0, 1], r=2,
                                       flatten=True, include_input=True)
        assert np.size(a) == 17
        assert np.all(np.in1d([0, 1], a))

    def test_is_fully_connected(self):
        assert self.net._is_fully_connected()
        # Disconnect pore 0 from the network
        Ps = self.net.find_neighbor_pores(pores=0)
        Ts = self.net.find_neighbor_throats(pores=0)
        op.topotools.trim(self.net, throats=Ts)
        assert not self.net._is_fully_connected()
        # Reconnect pore 0 to the network
        op.topotools.connect_pores(self.net, pores1=0, pores2=Ps)


if __name__ == '__main__':

    t = GenericNetworkTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
