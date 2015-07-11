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

    def test_connectivity_invalid(self):
        flag = False
        try:
            OpenPNM.Network.Cubic(shape=[3, 3, 3], connectivity=25)
        except:
            flag = True
        assert flag

    def test_shape_2D(self):
        net = OpenPNM.Network.Cubic(shape=[3, 3])
        assert net._shape == (3, 3, 1)

    def test_asarray(self):
        nums = self.net.Ps
        arr = self.net.asarray(nums)
        assert sp.shape(arr) == (3, 3, 3)

    def test_asarray_too_short(self):
        nums = self.net.Ps
        flag = False
        try:
            self.net.asarray(nums[:-1])
        except:
            flag = True
        assert flag

    def test_asarray_too_long(self):
        nums = sp.ones((self.net.Np+1,))
        flag = False
        try:
            self.net.asarray(nums)
        except:
            flag = True
        assert flag

    def test_from_array(self):
        arr = sp.ones(self.net._shape, dtype=int)
        self.net.fromarray(array=arr, propname='pore.test')
        assert sp.sum(self.net['pore.test']) == self.net.Np

    def test_from_array_wrong_size(self):
        arr = sp.ones(self.net._shape[:-1], dtype=int)
        flag = False
        try:
            self.net.fromarray(array=arr, propname='pore.test')
        except:
            flag = True
        assert flag

    def test_domain_area(self):
        A = self.net.domain_area(face=self.net.pores('top'))
        assert sp.allclose(A, 25, rtol=0.1)

    def test_domain_length(self):
        L = self.net.domain_length(face_1=self.net.pores('top'),
                                   face_2=self.net.pores('bottom'))
        assert sp.allclose(L, 4, rtol=1e-02)
