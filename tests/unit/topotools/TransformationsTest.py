import openpnm as op
import numpy as np
from numpy.testing import assert_allclose
from openpnm import topotools


class TransformationsTest:

    def setup_class(self):
        self.ws = op.Workspace()

    def teardown_class(self):
        self.ws.clear()

    def test_rotate_points(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        topotools.rotate_coords(network=net, a=45)
        a = np.array([[5.00000000e-01,  0.00000000e+00, 7.07106781e-01],
                      [5.00000000e-01, -7.07106781e-01, 1.41421356e+00]])
        assert np.allclose(net['pore.coords'][0:2, :], a)
        net = op.network.Cubic(shape=[3, 3, 3])
        topotools.rotate_coords(network=net, a=45, b=45, c=45)
        b = np.array([[-0.10355339, -0.10355339, 0.85355339],
                      [+0.04289322, -0.95710678, 1.35355339]])
        assert np.allclose(net['pore.coords'][0:2, :], b)
        net = op.network.Cubic(shape=[3, 3, 3])
        q = 45
        R = np.array([[1, 0, 0],
                      [0, np.cos(np.deg2rad(q)), -np.sin(np.deg2rad(q))],
                      [0, np.sin(np.deg2rad(q)), np.cos(np.deg2rad(q))]])
        topotools.rotate_coords(network=net, R=R)
        assert np.allclose(net['pore.coords'][0:2, :], a)

    def test_shear_points(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        a = np.array([[0.5, 0.5, 0.5],
                      [0.5, 0.5, 1.5]])
        assert np.allclose(net['pore.coords'][0:2, :], a)
        topotools.shear_coords(network=net, ay=1)
        b = np.array([[1.0, 0.5, 0.5],
                      [1.0, 0.5, 1.5]])
        assert np.allclose(net['pore.coords'][0:2, :], b)
        net = op.network.Cubic(shape=[3, 3, 3])
        S = [[1.0, 1.0, 0.0],
             [0.0, 1.0, 0.0],
             [0.0, 0.0, 1.0]]
        topotools.shear_coords(network=net, S=S)
        assert np.allclose(net['pore.coords'][0:2, :], b)


if __name__ == '__main__':

    t = TransformationsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
