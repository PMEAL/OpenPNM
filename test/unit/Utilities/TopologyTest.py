import OpenPNM
from OpenPNM.Utilities import topology
import scipy as sp
topo = topology()


class TopologyTest:

    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5], spacing=1)
        Ps = self.net.pores()
        Ts = self.net.throats()
        self.geom = OpenPNM.Geometry.TestGeometry(network=self.net,
                                                  pores=Ps,
                                                  throats=Ts,
                                                  name='test_geom')

    def test_extend(self):
        pass

    def test_trim_occluded_throats(self):
        Np = self.net.num_pores()
        Nt = self.net.num_throats()
        ts = self.net.find_neighbor_throats(pores=[0])
        self.geom['throat.area'][ts] = 0.0
        topo.trim_occluded_throats(self.net)
        assert self.net.num_pores() == Np - 1
        assert self.net.num_throats() == Nt - sp.size(ts)

    def test_merge_pores(self):
        net = OpenPNM.Network.Cubic(shape=[5, 5, 5], spacing=1)
        P = net.find_nearby_pores(pores=67, distance=2, flatten=True)
        topo.merge_pores(network=net, pores=P, labels=['merged'])
        assert net.Np == 95
        assert net.Nt == 222

if __name__ == '__main__':
    self = TopologyTest()
    self.setup_class()

