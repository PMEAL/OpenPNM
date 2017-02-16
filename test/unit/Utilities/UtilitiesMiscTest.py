import OpenPNM
import OpenPNM.Utilities.misc as misc
import scipy as sp
import time


class UtilitiesMiscTest:

    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[10, 10, 10])
        self.air = OpenPNM.Phases.Air(network=self.net)

    def test_find_path_single_pore_pair(self):
        a = misc.find_path(network=self.net, pore_pairs=[0, 1])
        assert sorted(list(a.keys())) == ['pores', 'throats']
        assert len(a['pores']) == 1
        assert len(a['throats']) == 1
        assert len(a['pores'][0]) == 2
        assert len(a['throats'][0]) == 1

    def test_find_path_multiple_pore_pairs(self):
        a = misc.find_path(network=self.net,
                           pore_pairs=([0, 1], [3, 6], [0, 8]))
        assert len(a['pores']) == 3
        assert len(a['throats']) == 3

    def test_find_path_with_weights(self):
        w = sp.ones_like(self.net.Ts)
        w[0] = self.net.Nt + 1
        a = misc.find_path(network=self.net, pore_pairs=([0, 1]), weights=w)
        assert len(a['pores'][0]) > 2
        assert len(a['throats'][0]) > 1

    def test_amalgamate_data(self):
        dict_ = misc.amalgamate_data(objs=[self.net, self.air])
        assert 'pore.'+self.air.name+'_molecular_weight' in dict_.keys()
        dict_ = misc.amalgamate_data(objs=[self.net, self.air], delimiter='|')
        assert 'pore.'+self.air.name+'|molecular_weight' in dict_.keys()

    def test_tic_toc(self):
        t = misc.tic()
        assert t is None
        t = misc.toc(quiet=False)
        assert t is None
        time.sleep(0.005)
        t = misc.toc(quiet=True)
        assert t > 0

    def test_iscoplanar(self):
        net = OpenPNM.Network.Cubic(shape=[2, 2, 2])
        pts1 = net.pores(labels=['top', 'left'], mode='intersection')
        pts2 = net.pores(labels=['bottom', 'right'], mode='intersection')
        pts = sp.hstack([pts1, pts2])
        check = misc.iscoplanar(coords=net['pore.coords'][pts])
        assert check
