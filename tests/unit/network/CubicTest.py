import openpnm as op
import scipy as sp
import pytest


class CubicTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_spacing_2D(self):
        net = op.network.Cubic(shape=[5, 5, 1], spacing=[1, 1])
        assert sp.all(net.spacing == [1.0, 1.0, 0.0])

    def test_spacing_3D(self):
        net = op.network.Cubic(shape=[5, 5, 5], spacing=[1, 1, 1])
        assert sp.all(net.spacing == [1.0, 1.0, 1.0])

    def test_spacing_2D_uneven(self):
        net = op.network.Cubic(shape=[5, 5, 1], spacing=[1, 2])
        assert sp.all(net.spacing == [1.0, 2.0, 0.0])

    def test_spacing_3D_uneven(self):
        net = op.network.Cubic(shape=[3, 4, 5], spacing=[1, 2, 3])
        assert sp.all(net.spacing == [1.0, 2.0, 3.0])

    def test_shape_1D(self):
        net = op.network.Cubic(shape=[5, 1, 1])
        assert sp.all(net.shape == [5, 1, 1])
        net = op.network.Cubic(shape=[1, 1, 5])
        assert sp.all(net.shape == [1, 1, 5])
        net = op.network.Cubic(shape=[1, 5, 1])
        assert sp.all(net.shape == [1, 5, 1])

    def test_shape_2D(self):
        net = op.network.Cubic(shape=[5, 5, 1])
        assert sp.all(net.shape == [5, 5, 1])
        net = op.network.Cubic(shape=[5, 1, 5])
        assert sp.all(net.shape == [5, 1, 5])
        net = op.network.Cubic(shape=[1, 5, 5])
        assert sp.all(net.shape == [1, 5, 5])

    def test_shape_3D(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        assert sp.all(net.shape == [5, 5, 5])
        net = op.network.Cubic(shape=[3, 4, 5])
        assert sp.all(net.shape == [3, 4, 5])
        net = op.network.Cubic(shape=[1, 5, 1])
        assert sp.all(net.shape == [1, 5, 1])

    def test_spacing_3D_rotated(self):
        net = op.network.Cubic(shape=[5, 5, 5], spacing=[1, 1, 1])
        theta = 0.1
        R = sp.array([[1, 0, 0],
                      [0, sp.cos(theta), -sp.sin(theta)],
                      [0, sp.sin(theta), sp.cos(theta)]])
        net['pore.coords'] = sp.tensordot(net['pore.coords'], R, axes=(1, 1))
        assert sp.all(net.spacing == [1.0, 1.0, 1.0])

    def test_spacing_3D_rotated_uneven(self):
        net = op.network.Cubic(shape=[3, 4, 5], spacing=[1, 2, 3])
        theta = 0.1
        R = sp.array([[1, 0, 0],
                      [0, sp.cos(theta), -sp.sin(theta)],
                      [0, sp.sin(theta), sp.cos(theta)]])
        net['pore.coords'] = sp.tensordot(net['pore.coords'], R, axes=(1, 1))
        assert sp.all(net.spacing == [1.0, 2.0, 3.0])

    def test_spacing_2D_sheared(self):
        net = op.network.Cubic(shape=[5, 5, 1], spacing=1)
        S = sp.array([[1, 1, 0],
                      [0, 1, 0],
                      [0, 0, 1]])
        net['pore.coords'] = (S@net['pore.coords'].T).T
        assert sp.allclose(net.spacing, [1.0, 2**0.5, 0.0])

    def test_spacing_2D_sheared_uneven(self):
        net = op.network.Cubic(shape=[5, 5, 1], spacing=[1, 2])
        S = sp.array([[1, 1, 0],
                      [0, 1, 0],
                      [0, 0, 1]])
        net['pore.coords'] = (S@net['pore.coords'].T).T
        assert sp.allclose(net.spacing, [1.0, 2*(2**0.5), 0.0])

    def test_spacing_3D_sheared(self):
        net = op.network.Cubic(shape=[5, 5, 3], spacing=1)
        S = sp.array([[1, 1, 0],
                      [0, 1, 0],
                      [0, 0, 1]])
        net['pore.coords'] = (S@net['pore.coords'].T).T
        assert sp.allclose(net.spacing, [1.0, 2**0.5, 1.0])

    def test_spacing_3D_sheared_uneven(self):
        net = op.network.Cubic(shape=[3, 4, 5], spacing=[1, 2, 3])
        S = sp.array([[1, 1, 0],
                      [0, 1, 0],
                      [0, 0, 1]])
        net['pore.coords'] = (S@net['pore.coords'].T).T
        assert sp.allclose(net.spacing, [1.0, 2*(2**0.5), 3.0])

    def test_spacing_on_joggled_network(self):
        net = op.network.Cubic(shape=[3, 4, 5])
        net['pore.coords'] += sp.rand(net.Np, 3)
        with pytest.raises(Exception):
            net.spacing

    def test_spacing_on_network_with_boundary_pores(self):
        net = op.network.Cubic(shape=[3, 4, 5])
        net.add_boundary_pores()
        with pytest.raises(Exception):
            net.spacing


if __name__ == '__main__':

    t = CubicTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
