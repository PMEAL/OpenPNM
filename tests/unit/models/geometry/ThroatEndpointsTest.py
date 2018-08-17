import openpnm as op
import openpnm.models.geometry.throat_endpoints as mods
import numpy as np
from numpy.testing import assert_allclose


class ThroatEndpointsTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[1, 2, 1], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.base = np.array([0.5, 0.5, 0.5])

    def test_spherical_pores(self):
        self.geo['pore.diameter'] = 0.5
        self.geo['throat.diameter'] = 0.25
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.spherical_pores,
                           regen_mode='normal')
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.2165064, 0]) + self.base
        EP2d = np.array([0, 1 - 0.2165064, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_cubic_pores(self):
        self.geo['pore.diameter'] = 0.5
        self.geo['throat.diameter'] = 0.25
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.cubic_pores,
                           regen_mode='normal')
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.25, 0]) + self.base
        EP2d = np.array([0, 1 - 0.25, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_spherical_pores_with_apparent_overlap(self):
        # Apparent overlap means pores are overlapping, but since throat
        # diameter is large enough to encompass the pores' intersection area
        self.geo['pore.diameter'] = np.array([1.3, 1.5])
        self.geo['throat.diameter'] = 1.1
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.spherical_pores,
                           regen_mode='normal')
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.3464102, 0]) + self.base
        EP2d = np.array([0, 1 - 0.5099020, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_spherical_pores_with_true_overlap(self):
        self.geo['pore.diameter'] = np.array([1.3, 1.5])
        self.geo['throat.diameter'] = 0.6
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.spherical_pores,
                           regen_mode='normal')
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.43, 0]) + self.base
        EP2d = np.array([0, 0.43, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_cubic_pores_with_overlap(self):
        self.geo['pore.diameter'] = np.array([1.3, 1.5])
        self.geo['throat.diameter'] = 1.1
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.cubic_pores,
                           regen_mode='normal')
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.25, 0]) + self.base
        EP2d = np.array([0, 0.25, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_circular_pores(self):
        self.geo['pore.diameter'] = 0.5
        self.geo['throat.diameter'] = 0.25
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.spherical_pores,
                           regen_mode='normal')
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.2165064, 0]) + self.base
        EP2d = np.array([0, 1 - 0.2165064, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_circular_pores_with_apparent_overlap(self):
        # Apparent overlap means pores are overlapping, but since throat
        # diameter is large enough to encompass the pores' intersection area
        self.geo['pore.diameter'] = np.array([1.3, 1.5])
        self.geo['throat.diameter'] = 1.1
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.circular_pores,
                           regen_mode='normal')
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.3464102, 0]) + self.base
        EP2d = np.array([0, 1 - 0.5099020, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_circular_pores_with_true_overlap(self):
        self.geo['pore.diameter'] = np.array([1.3, 1.5])
        self.geo['throat.diameter'] = 0.6
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.circular_pores,
                           regen_mode='normal')
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.43, 0]) + self.base
        EP2d = np.array([0, 0.43, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_spherical_pores_with_zero_diameter(self):
        self.geo['pore.diameter'] = np.array([0.5, 0.0])
        self.geo['throat.diameter'] = 0.25
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.spherical_pores,
                           regen_mode='normal')
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.2165064, 0]) + self.base
        EP2d = np.array([0, 1, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_circular_pores_with_zero_diameter(self):
        self.geo['pore.diameter'] = np.array([0.5, 0.0])
        self.geo['throat.diameter'] = 0.25
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.circular_pores,
                           regen_mode='normal')
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.2165064, 0]) + self.base
        EP2d = np.array([0, 1, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_cubic_pores_with_zero_diameter(self):
        self.geo['pore.diameter'] = np.array([0.5, 0.0])
        self.geo['throat.diameter'] = 0.25
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.cubic_pores,
                           regen_mode='normal')
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.25, 0]) + self.base
        EP2d = np.array([0, 1, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

if __name__ == '__main__':

    t = ThroatEndpointsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
