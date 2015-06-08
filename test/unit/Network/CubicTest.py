import OpenPNM
import scipy as sp


class CubicTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5], spacing=1)
        self.net['pore.diameter'] = sp.rand(self.net.Np)

    def test_connectivity_6(self):
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3], connectivity=6)
        assert net.num_neighbors(pores=13) == 6

    def test_connectivity_8(self):
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3], connectivity=8)
        assert net.num_neighbors(pores=13) == 8

    def test_connectivity_12(self):
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3], connectivity=12)
        assert net.num_neighbors(pores=13) == 12

    def test_connectivity_14(self):
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3], connectivity=14)
        assert net.num_neighbors(pores=13) == 14

    def test_connectivity_18(self):
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3], connectivity=18)
        assert net.num_neighbors(pores=13) == 18

    def test_connectivity_20(self):
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3], connectivity=20)
        assert net.num_neighbors(pores=13) == 20

    def test_connectivity_26(self):
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3], connectivity=26)
        assert net.num_neighbors(pores=13) == 26

    def test_domain_area(self):
        A = self.net.domain_area(face=self.net.pores('top'))
        assert sp.allclose(A, 25, rtol=1e-02)

    def test_domain_length(self):
        L = self.net.domain_length(face_1=self.net.pores('top'),
                                   face_2=self.net.pores('bottom'))
        assert sp.allclose(L, 4, rtol=1e-02)

    def test_from_array(self):
        arr = sp.ones(self.net._shape, dtype=int)
        self.net.fromarray(array=arr, propname='pore.test')
        assert sp.sum(self.net['pore.test']) == self.net.Np
