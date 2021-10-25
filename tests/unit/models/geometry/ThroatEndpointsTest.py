import openpnm as op
import openpnm.models.geometry.throat_endpoints as mods
import openpnm.models.geometry as gm
import openpnm.models.physics as pm
import numpy as np
from numpy.testing import assert_allclose


class ThroatEndpointsTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[1, 3, 1], spacing=1.0)
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
        # Throat 1->2
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.2165064, 0]) + self.base
        EP2d = np.array([0, 1 - 0.2165064, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)
        # Throat 2->3
        EP1 = self.geo['throat.endpoints.head'][1]
        EP2 = self.geo['throat.endpoints.tail'][1]
        EP1d = np.array([0, 1 + 0.2165064, 0]) + self.base
        EP2d = np.array([0, 2 - 0.2165064, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_spherical_pores_with_throat_centroid(self):
        self.geo['throat.centroid'] = np.array([[0, 0.5, 0.5],
                                               [0, 1.5, 0]]) + self.base
        self.geo['pore.diameter'] = 0.5
        self.geo['throat.diameter'] = 0.25
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.spherical_pores,
                           regen_mode='normal')
        # Throat 1->2
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.2165064 * np.cos(np.pi/4),
                         0.2165064 * np.sin(np.pi/4)]) + self.base
        EP2d = np.array([0, 1 - 0.2165064 * np.cos(np.pi/4),
                         0.2165064 * np.sin(np.pi/4)]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)
        # Throat 2->3
        EP1 = self.geo['throat.endpoints.head'][1]
        EP2 = self.geo['throat.endpoints.tail'][1]
        EP1d = np.array([0, 1 + 0.2165064, 0]) + self.base
        EP2d = np.array([0, 2 - 0.2165064, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)
        del self.geo['throat.centroid']

    def test_spherical_pores_with_throat_centroid_and_overlap(self):
        self.geo['throat.centroid'] = np.array([[0, 0.5, 1],
                                                [0, 1.5, 0]]) + self.base
        self.geo['pore.diameter'] = [1.5, 1.0, 0.5]
        self.geo['throat.diameter'] = 0.25
        self.geo.add_model(propname='throat.endpoints',
                            model=mods.spherical_pores,
                            regen_mode='normal')
        # Only check throat 1->2, 2->3 is already tested (no-overlap)
        EP12_head = self.geo['throat.endpoints.head'][0]
        EP12_tail = self.geo['throat.endpoints.tail'][0]
        EP12_head_desired = np.array([0, 0.29348392, 0.58696784]) + self.base
        EP12_tail_desired = np.array([0, 0.84627033, 0.30745935]) + self.base
        assert_allclose(EP12_head, desired=EP12_head_desired)
        assert_allclose(EP12_tail, desired=EP12_tail_desired)
        del self.geo["throat.centroid"]

    def test_cubic_pores(self):
        self.geo['pore.diameter'] = 0.5
        self.geo['throat.diameter'] = 0.25
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.cubic_pores,
                           regen_mode='normal')
        # Throat 1->2
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.25, 0]) + self.base
        EP2d = np.array([0, 1 - 0.25, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)
        # Throat 2->3
        EP1 = self.geo['throat.endpoints.head'][1]
        EP2 = self.geo['throat.endpoints.tail'][1]
        EP1d = np.array([0, 1 + 0.25, 0]) + self.base
        EP2d = np.array([0, 2 - 0.25, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_square_pores(self):
        self.geo['pore.diameter'] = 0.5
        self.geo['throat.diameter'] = 0.25
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.square_pores,
                           regen_mode='normal')
        # Throat 1->2
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.25, 0]) + self.base
        EP2d = np.array([0, 1 - 0.25, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)
        # Throat 2->3
        EP1 = self.geo['throat.endpoints.head'][1]
        EP2 = self.geo['throat.endpoints.tail'][1]
        EP1d = np.array([0, 1 + 0.25, 0]) + self.base
        EP2d = np.array([0, 2 - 0.25, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_spherical_pores_with_apparent_overlap(self):
        # Apparent overlap means pores are overlapping, but since throat
        # diameter is large enough, it encompasses the pores' intersection area
        self.geo['pore.diameter'] = np.array([1.3, 1.5, 1.3])
        self.geo['throat.diameter'] = 1.1
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.spherical_pores,
                           regen_mode='normal')
        # Throat 1->2
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.3464102, 0]) + self.base
        EP2d = np.array([0, 1 - 0.5099020, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)
        # Throat 2->3
        EP1 = self.geo['throat.endpoints.head'][1]
        EP2 = self.geo['throat.endpoints.tail'][1]
        EP1d = np.array([0, 1 + 0.5099020, 0]) + self.base
        EP2d = np.array([0, 2 - 0.3464102, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_spherical_pores_with_true_overlap(self):
        self.geo['pore.diameter'] = np.array([1.3, 1.5, 1.3])
        self.geo['throat.diameter'] = 0.6
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.spherical_pores,
                           regen_mode='normal')
        # Throat 1->2
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.43, 0]) + self.base
        EP2d = np.array([0, 0.43, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)
        # Throat 2->3
        EP1 = self.geo['throat.endpoints.head'][1]
        EP2 = self.geo['throat.endpoints.tail'][1]
        EP1d = np.array([0, 1 + 0.57, 0]) + self.base
        EP2d = np.array([0, 1 + 0.57, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_cubic_pores_with_overlap(self):
        self.geo['pore.diameter'] = np.array([1.3, 1.5, 1.3])
        self.geo['throat.diameter'] = 1.1
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.cubic_pores,
                           regen_mode='normal')
        # Throat 1->2
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.25, 0]) + self.base
        EP2d = np.array([0, 0.25, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)
        # Throat 2->3
        EP1 = self.geo['throat.endpoints.head'][1]
        EP2 = self.geo['throat.endpoints.tail'][1]
        EP1d = np.array([0, 1 + 0.75, 0]) + self.base
        EP2d = np.array([0, 1 + 0.75, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_circular_pores(self):
        self.geo['pore.diameter'] = 0.5
        self.geo['throat.diameter'] = 0.25
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.spherical_pores,
                           regen_mode='normal')
        # Throat 1->2
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.2165064, 0]) + self.base
        EP2d = np.array([0, 1 - 0.2165064, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)
        # Throat 2->3
        EP1 = self.geo['throat.endpoints.head'][1]
        EP2 = self.geo['throat.endpoints.tail'][1]
        EP1d = np.array([0, 1 + 0.2165064, 0]) + self.base
        EP2d = np.array([0, 2 - 0.2165064, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_circular_pores_with_apparent_overlap(self):
        # Apparent overlap means pores are overlapping, but since throat
        # diameter is large enough, it encompasses the pores' intersection area
        self.geo['pore.diameter'] = np.array([1.3, 1.5, 1.3])
        self.geo['throat.diameter'] = 1.1
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.circular_pores,
                           regen_mode='normal')
        # Throat 1->2
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.3464102, 0]) + self.base
        EP2d = np.array([0, 1 - 0.5099020, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)
        # Throat 2->3
        EP1 = self.geo['throat.endpoints.head'][1]
        EP2 = self.geo['throat.endpoints.tail'][1]
        EP1d = np.array([0, 1 + 0.5099020, 0]) + self.base
        EP2d = np.array([0, 2 - 0.3464102, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_circular_pores_with_true_overlap(self):
        self.geo['pore.diameter'] = np.array([1.3, 1.5, 1.3])
        self.geo['throat.diameter'] = 0.6
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.circular_pores,
                           regen_mode='normal')
        # Throat 1->2
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.43, 0]) + self.base
        EP2d = np.array([0, 0.43, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)
        # Throat 2->3
        EP1 = self.geo['throat.endpoints.head'][1]
        EP2 = self.geo['throat.endpoints.tail'][1]
        EP1d = np.array([0, 1 + 0.57, 0]) + self.base
        EP2d = np.array([0, 1 + 0.57, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_spherical_pores_with_zero_diameter(self):
        self.geo['pore.diameter'] = np.array([0.5, 0.0, 0.5])
        self.geo['throat.diameter'] = 0.25
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.spherical_pores,
                           regen_mode='normal')
        # Throat 1->2
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.2165064, 0]) + self.base
        EP2d = np.array([0, 1, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)
        # Throat 2->3
        EP1 = self.geo['throat.endpoints.head'][1]
        EP2 = self.geo['throat.endpoints.tail'][1]
        EP1d = np.array([0, 1, 0]) + self.base
        EP2d = np.array([0, 2 - 0.2165064, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_circular_pores_with_zero_diameter(self):
        self.geo['pore.diameter'] = np.array([0.5, 0.0, 0.5])
        self.geo['throat.diameter'] = 0.25
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.circular_pores,
                           regen_mode='normal')
        # Throat 1->2
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.2165064, 0]) + self.base
        EP2d = np.array([0, 1, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)
        # Throat 2->3
        EP1 = self.geo['throat.endpoints.head'][1]
        EP2 = self.geo['throat.endpoints.tail'][1]
        EP1d = np.array([0, 1, 0]) + self.base
        EP2d = np.array([0, 2 - 0.2165064, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)

    def test_cubic_pores_with_zero_diameter(self):
        self.geo['pore.diameter'] = np.array([0.5, 0.0, 0.5])
        self.geo['throat.diameter'] = 0.25
        self.geo.add_model(propname='throat.endpoints',
                           model=mods.cubic_pores,
                           regen_mode='normal')
        # Throat 1->2
        EP1 = self.geo['throat.endpoints.head'][0]
        EP2 = self.geo['throat.endpoints.tail'][0]
        EP1d = np.array([0, 0.25, 0]) + self.base
        EP2d = np.array([0, 1, 0]) + self.base
        assert_allclose(EP1, desired=EP1d)
        assert_allclose(EP2, desired=EP2d)
        # Throat 2->3
        EP1 = self.geo['throat.endpoints.head'][1]
        EP2 = self.geo['throat.endpoints.tail'][1]
        EP1d = np.array([0, 1, 0]) + self.base
        EP2d = np.array([0, 2 - 0.25, 0]) + self.base
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
