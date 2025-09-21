import numpy as np
import pytest

import openpnm as op
from openpnm.topotools import extend, trim


class CheckNetworkHealthTest:

    def setup_class(self):
        pass

    def test_check_network_health_healthy(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        a = op.utils.check_network_health(net)
        items = set(['headless_throats',
                     'looped_throats',
                     'isolated_pores',
                     'disconnected_pores',
                     'duplicate_throats',
                     'bidirectional_throats'])
        assert items == a.keys()
        assert np.size(list(a.values())) == 0

    def test_check_network_health_one_isolated_cluster(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        Ps = net['pore.coords'][:, 2] == 2.5
        Ps = Ps*(net['pore.coords'][:, 1] > 2.5)
        Ts = net.find_neighbor_throats(pores=Ps, mode='exclusive_or')
        trim(network=net, throats=Ts)
        a = op.utils.check_network_health(net)
        assert len(a['disconnected_pores']) == 10

    def test_check_network_health_two_isolated_clusters(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        # Create first isolated cluster
        Ps = net['pore.coords'][:, 2] == 2.5
        Ps = Ps*(net['pore.coords'][:, 1] > 2.5)
        Ts = net.find_neighbor_throats(pores=Ps, mode='exclusive_or')
        trim(network=net, throats=Ts)
        # Create a second isolated cluster
        Ps = net['pore.coords'][:, 2] == 3.5
        Ps = Ps*(net['pore.coords'][:, 1] < 2.5)
        Ts = net.find_neighbor_throats(pores=Ps, mode='exclusive_or')
        trim(network=net, throats=Ts)
        a = op.utils.check_network_health(net)
        # assert len(a['disconnected_clusters']) == 3
        assert len(a['disconnected_pores']) == 20

    def test_check_network_health_isolated_pores_and_clusters(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        # Create first isolated cluster
        Ps = net['pore.coords'][:, 2] == 2.5
        Ps = Ps*(net['pore.coords'][:, 1] > 2.5)
        Ts = net.find_neighbor_throats(pores=Ps, mode='exclusive_or')
        trim(network=net, throats=Ts)
        # Create an isolated pore
        Ts = net.find_neighbor_throats(pores=0)
        trim(network=net, throats=Ts)
        a = op.utils.check_network_health(net)
        # Ensure trim_pores has right length
        assert len(a['disconnected_pores']) == 11
        # Ensure 0 is listed in trim pores
        assert 0 in a['disconnected_pores']

    def test_check_network_health_isolated_pores(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        Ts = net.find_neighbor_throats(pores=0)
        trim(network=net, throats=Ts)
        a = op.utils.check_network_health(net)
        assert a['isolated_pores'] == np.array([0])
        trim(network=net, pores=a['disconnected_pores'])
        a = op.utils.check_network_health(net)
        assert np.size(a['disconnected_pores']) == 0

    def test_check_network_health_duplicate_throat(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        P12 = net['throat.conns'][0]
        extend(network=net, throat_conns=[P12])
        a = op.utils.check_network_health(net)
        assert len(a['duplicate_throats']) == 1
        assert a['duplicate_throats'][0] == 300

    def test_check_network_health_triplicate_throats(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        P12 = net['throat.conns'][0]
        extend(network=net, throat_conns=[P12])
        extend(network=net, throat_conns=[P12])
        a = op.utils.check_network_health(net)
        assert len(a['duplicate_throats']) == 2

    def test_check_network_health_multiple_duplicate_throats(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        P12 = net['throat.conns'][0]
        extend(network=net, throat_conns=[P12])
        P12 = net['throat.conns'][1]
        extend(network=net, throat_conns=[P12])
        a = op.utils.check_network_health(net)
        assert len(a['duplicate_throats']) == 2

    def test_check_network_health_bidirectional_throats(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        P12 = net['throat.conns'][0]
        net['throat.conns'][0] = [P12[1], P12[0]]
        a = op.utils.check_network_health(net)
        assert np.size(a['bidirectional_throats']) == 1
        assert np.size(a['duplicate_throats']) == 0

    # def test_check_network_health_headless_throats(self):
    #     net = op.network.Cubic(shape=[5, 5, 5])
    #     with pytest.raises(Exception):
    #         extend(network=net, throat_conns=[[5, 5555]])
    #     net['throat.conns'][0] = [5, 5555]
    #     a = net.check_network_health()
    #     assert a['headless_throats'] == np.array([0])

    def test_check_network_health_looped_throats(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        extend(network=net, throat_conns=[[5, 5]])
        a = op.utils.check_network_health(net)
        assert a['looped_throats'] == np.array([300])


if __name__ == '__main__':

    t = CheckNetworkHealthTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            t.__getattribute__(item)()
