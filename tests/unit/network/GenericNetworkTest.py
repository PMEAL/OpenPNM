import openpnm as op
import scipy as sp


class GenericNetworkTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[10, 10, 10])

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_find_connected_pores_numeric_not_flattend(self):
        a = self.net.find_connected_pores(throats=[0, 1])
        assert sp.all(a.flatten() == [0, 1, 1, 2])

    def test_find_connected_pores_numeric_flattend(self):
        a = self.net.find_connected_pores(throats=[0, 1], flatten=True)
        assert sp.all(a == [0, 1, 2])

    def test_find_connected_pores_boolean_flattend(self):
        Tind = sp.zeros((self.net.Nt,), dtype=bool)
        Tind[[0, 1]] = True
        a = self.net.find_connected_pores(throats=Tind, flatten=True)
        assert sp.all(a == [0, 1, 2])

    def test_find_connected_pores_empty_flattend(self):
        a = self.net.find_connected_pores(throats=[], flatten=True)
        assert sp.shape(a) == (0, )

    def test_find_neighbor_pores_numeric(self):
        a = self.net.find_neighbor_pores(pores=[])
        assert sp.size(a) == 0

    def test_find_neighbor_pores_boolean(self):
        Pind = sp.zeros((self.net.Np,), dtype=bool)
        Pind[[0, 1]] = True
        a = self.net.find_neighbor_pores(pores=Pind)
        assert sp.all(a == [2, 10, 11, 100, 101])

    def test_find_neighbor_pores_numeric_union(self):
        a = self.net.find_neighbor_pores(pores=[0, 2], mode='union')
        assert sp.all(a == [1, 3, 10, 12, 100, 102])

    def test_find_neighbor_pores_numeric_intersection(self):
        a = self.net.find_neighbor_pores(pores=[0, 2], mode='xnor')
        assert sp.all(a == [1])

    def test_find_neighbor_pores_numeric_exclusive_or(self):
        a = self.net.find_neighbor_pores(pores=[0, 2], mode='exclusive_or')
        assert sp.all(a == [3, 10, 12, 100, 102])

    def test_find_neighbor_pores_numeric_union_incl_self(self):
        a = self.net.find_neighbor_pores(pores=[0, 2], mode='or',
                                         exclude_input=False)
        assert sp.all(a == [0, 1, 2, 3, 10, 12, 100, 102])

    def test_find_neighbor_pores_numeric_intersection_incl_self(self):
        a = self.net.find_neighbor_pores(pores=[0, 2], mode='and',
                                         exclude_input=False)
        assert sp.all(a == [0, 1, 2])

    def test_find_neighbor_pores_numeric_intersection_excl_self(self):
        a = self.net.find_neighbor_pores(pores=[0, 2], mode='and',
                                         exclude_input=True)
        assert sp.all(a == [1])

    def test_find_neighbor_pores_numeric_exclusive_or_incl_self(self):
        a = self.net.find_neighbor_pores(pores=[0, 2],
                                         mode='exclusive_or',
                                         exclude_input=False)
        assert sp.all(a == [0, 2, 3, 10, 12, 100, 102])

    def test_find_neighbor_throats_empty(self):
        a = self.net.find_neighbor_throats(pores=[])
        assert sp.size(a) == 0

    def test_find_neighbor_throats_boolean(self):
        Pind = sp.zeros((self.net.Np,), dtype=bool)
        Pind[[0, 1]] = True
        a = self.net.find_neighbor_throats(pores=Pind)
        assert sp.all(a == [0, 1, 900, 901, 1800, 1801])

    def test_find_neighbor_throats_numeric_union(self):
        a = self.net.find_neighbor_throats(pores=[0, 2], mode='union')
        assert sp.all(a == [0, 1, 2, 900, 902, 1800, 1802])

    def test_find_neighbor_throats_numeric_intersection(self):
        a = self.net.find_neighbor_throats(pores=[0, 2], mode='xnor')
        assert sp.size(a) == 0

    def test_find_neighbor_throats_numeric_exclusive_or(self):
        a = self.net.find_neighbor_throats(pores=[0, 2],
                                           mode='exclusive_or')
        assert sp.all(a == [0, 1, 2, 900, 902, 1800, 1802])

    def test_num_neighbors_empty(self):
        a = self.net.num_neighbors(pores=[])
        assert sp.size(a) == 0

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
        assert sp.all(a == [3, 4])
        a = self.net.num_neighbors(pores=0, flatten=False)
        assert sp.all(a == [3])
        assert isinstance(a, sp.ndarray)

    def test_find_nearby_pores_distance_1(self):
        a = self.net.find_nearby_pores(pores=[0, 1], r=1, exclude_input=True)
        b = self.net.find_neighbor_pores(pores=[0, 1], flatten=False,
                                         exclude_input=False)
        assert sp.all([sp.all(a[i] == b[i]) for i in range(0, len(a))])

    def test_find_nearby_pores_distance_2(self):
        a = self.net.find_nearby_pores(pores=[0, 1], r=2)
        assert sp.all([sp.size(a[i]) for i in [0, 1]] == [10, 14])

    def test_find_nearby_pores_distance_0(self):
        a = self.net.find_nearby_pores(pores=[0, 1], r=1e-9)
        assert sp.shape(a) == (2, 0)

    def test_find_nearby_pores_distance_1_flattened(self):
        a = self.net.find_nearby_pores(pores=[0, 1], r=1, flatten=True)
        b = self.net.find_neighbor_pores(pores=[0, 1])
        assert sp.all(a == b)

    def test_find_nearby_pores_distance_2_flattened(self):
        a = self.net.find_nearby_pores(pores=[0, 1], r=2, flatten=True)
        assert sp.size(a) == 15

    def test_find_nearby_pores_distance_2_flattened_inclself(self):
        a = self.net.find_nearby_pores(pores=[0, 1], r=2,
                                       flatten=True, exclude_input=False)
        assert sp.size(a) == 17
        assert sp.all(sp.in1d([0, 1], a))


if __name__ == '__main__':

    t = GenericNetworkTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
