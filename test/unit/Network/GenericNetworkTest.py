import OpenPNM
import scipy as sp


class GenericNetworkTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[10, 10, 10])

    def teardown_class(self):
        ctrl = OpenPNM.Base.Controller()
        ctrl.clear()

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
        assert sp.shape(a) == (0, 2)

    def test_find_neighbor_pores_numeric(self):
        a = self.net.find_neighbor_pores(pores=[])
        assert sp.size(a) == 0

    def test_find_neighbor_pores_boolean(self):
        Pind = sp.zeros((self.net.Np,), dtype=bool)
        Pind[[0, 1]] = True
        a = self.net.find_neighbor_pores(pores=Pind)
        assert sp.all(a == [2, 10, 11, 100, 101])

    def test_find_neighbor_pores_numeric_union(self):
        a = self.net.find_neighbor_pores(pores=[0, 2],
                                         mode='union')
        assert sp.all(a == [1, 3, 10, 12, 100, 102])

    def test_find_neighbor_pores_numeric_intersection(self):
        a = self.net.find_neighbor_pores(pores=[0, 2],
                                         mode='intersection')
        assert sp.all(a == [1])

    def test_find_neighbor_pores_numeric_notintersection(self):
        a = self.net.find_neighbor_pores(pores=[0, 2],
                                         mode='not_intersection')
        assert sp.all(a == [3, 10, 12, 100, 102])

    def test_find_neighbor_pores_numeric_union_incl_self(self):
        a = self.net.find_neighbor_pores(pores=[0, 2],
                                         mode='union',
                                         excl_self=False)
        assert sp.all(a == [0, 1, 2, 3, 10, 12, 100, 102])

    def test_find_neighbor_pores_numeric_intersection_incl_self(self):
        a = self.net.find_neighbor_pores(pores=[0, 2],
                                         mode='intersection',
                                         excl_self=False)
        assert sp.all(a == [1])

    def test_find_neighbor_pores_numeric_notintersection_incl_self(self):
        a = self.net.find_neighbor_pores(pores=[0, 2],
                                         mode='not_intersection',
                                         excl_self=False)
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
        a = self.net.find_neighbor_throats(pores=[0, 2], mode='intersection')
        assert sp.size(a) == 0

    def test_find_neighbor_throats_numeric_not_intersection(self):
        a = self.net.find_neighbor_throats(pores=[0, 2],
                                           mode='not_intersection')
        assert sp.all(a == [0, 1, 2, 900, 902, 1800, 1802])

    def test_num_neighbors_empty(self):
        a = self.net.num_neighbors(pores=[])
        assert sp.size(a) == 0

    def test_num_neighbors_boolean(self):
        Pind = sp.zeros((self.net.Np,), dtype=bool)
        Pind[0] = True
        a = self.net.num_neighbors(pores=Pind)
        assert a == 3

    def test_num_neighbors_numeric_flattened(self):
        a = self.net.num_neighbors(pores=[0, 2], flatten=True)
        assert a == 6
        assert isinstance(a, int)

    def test_num_neighbors_numeric_notflattened(self):
        a = self.net.num_neighbors(pores=[0, 2], flatten=False)
        assert sp.all(a == [3, 4])

    def test_num_neighbors_single_pore_notflattened(self):
        a = self.net.num_neighbors(pores=0, flatten=False)
        assert sp.all(a == [3])
        assert isinstance(a, sp.ndarray)

    def test_find_interface_throats(self):
        self.net['pore.domain1'] = False
        self.net['pore.domain2'] = False
        self.net['pore.domain3'] = False
        self.net['pore.domain1'][[0, 1, 2]] = True
        self.net['pore.domain2'][[100, 101, 102]] = True
        self.net['pore.domain3'][900:] = True
        a = self.net.find_interface_throats(labels=['domain1', 'domain2'])
        assert sp.all(a == [1800, 1801, 1802])
        a = self.net.find_interface_throats(labels=['domain1', 'domain3'])
        assert sp.size(a) == 0

    def test_check_network_health_healthy(self):
        a = self.net.check_network_health()
        items = set(['disconnected_clusters',
                     'isolated_pores',
                     'trim_pores',
                     'duplicate_throats',
                     'bidirectional_throats',
                     'headless_throats',
                     'looped_throats'])
        assert items == a.keys()
        assert sp.size(list(a.values())) == 0

    def test_check_network_isolated_pores(self):
        net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        Ts = net.find_neighbor_throats(pores=0)
        net.trim(throats=Ts)
        a = net.check_network_health()
        assert a['isolated_pores'] == 0
        net.trim(a['trim_pores'])
        a = net.check_network_health()
        assert sp.size(list(a.values())) == 0

    def test_check_network_health_duplicate_throat(self):
        net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        P12 = net['throat.conns'][0]
        net.extend(throat_conns=[P12])
        a = net.check_network_health()
        assert len(a['duplicate_throats'][0]) == 2

    def test_check_network_health_triplicate_throats(self):
        net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        P12 = net['throat.conns'][0]
        net.extend(throat_conns=[P12])
        net.extend(throat_conns=[P12])
        a = net.check_network_health()
        assert len(a['duplicate_throats'][0]) == 3

    def test_check_network_health_multiple_duplicate_throats(self):
        net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        P12 = net['throat.conns'][0]
        net.extend(throat_conns=[P12])
        P12 = net['throat.conns'][1]
        net.extend(throat_conns=[P12])
        a = net.check_network_health()
        assert len(a['duplicate_throats'][1]) == 2

    def test_check_network_health_bidirectional_throats(self):
        net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        P12 = net['throat.conns'][0]
        net['throat.conns'][0] = [P12[1], P12[0]]
        a = net.check_network_health()
        assert sp.size(a['bidirectional_throats']) == 1
        assert sp.size(a['duplicate_throats']) == 0

    def test_check_network_health_headless_throats(self):
        net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        net.extend(throat_conns=[[5, 5555]])
        a = net.check_network_health()
        assert a['headless_throats'] == sp.array([300])

    def test_check_network_health_looped_throats(self):
        net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        net.extend(throat_conns=[[5, 5]])
        a = net.check_network_health()
        assert a['looped_throats'] == sp.array([300])

    def test_find_nearby_pores_distance_1(self):
        a = self.net.find_nearby_pores(pores=[0, 1], distance=1)
        b = self.net.find_neighbor_pores(pores=[0, 1], flatten=False)
        assert sp.all([sp.all(a[i] == b[i]) for i in range(0, len(a))])

    def test_find_nearby_pores_distance_2(self):
        a = self.net.find_nearby_pores(pores=[0, 1], distance=2)
        assert sp.all([sp.size(a[i]) for i in [0, 1]] == [10, 14])

    def test_find_nearby_pores_distance_0(self):
        a = self.net.find_nearby_pores(pores=[0, 1], distance=0)
        assert sp.shape(a) == (2, 0)

    def test_find_nearby_pores_distance_1_flattened(self):
        a = self.net.find_nearby_pores(pores=[0, 1], distance=1, flatten=True)
        b = self.net.find_neighbor_pores(pores=[0, 1])
        assert sp.all(a == b)

    def test_find_nearby_pores_distance_2_flattened(self):
        a = self.net.find_nearby_pores(pores=[0, 1], distance=2, flatten=True)
        assert sp.size(a) == 15

    def test_find_nearby_pores_distance_2_flattened_inclself(self):
        a = self.net.find_nearby_pores(pores=[0, 1],
                                       distance=2,
                                       flatten=True,
                                       excl_self=False)
        assert sp.size(a) == 17
        assert sp.all(sp.in1d([0, 1], a))
