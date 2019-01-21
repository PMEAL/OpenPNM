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

    def test_straight_throat(self):
        r'''
        For a given network length and fixed cross-sectional area, check that
        the conductance is independent of number of pores and throat lengths
        '''
        # Domain length
        length = 5.0
        # Cross-sectional area of pores and throats (cylinder)
        # Results should be independent of this
        rad = 0.5
        area = np.pi*rad**2
        # Diffusivity in open air
        D_ab = 1.0
        # Number of internal pores in network
        N = np.arange(10, 110, 10)
        nets = []
        # Portion of the pore center to center distance attributed to throat
        # Any value below 1.0 should work
        throat_portion = 0.1

        def throat_centroid(network):
            cn = network['throat.conns']
            return np.mean(network['pore.coords'][cn], axis=1)

        def throat_vector(network):
            cn = network['throat.conns']
            vec = (network['pore.coords'][cn[:, 1]] -
                   network['pore.coords'][cn[:, 0]])
            return vec/np.linalg.norm(vec, axis=1)[:, np.newaxis]

        conds = []
        for n in N:
            lc = length/n
            net = op.network.Cubic(shape=[n, 1, 1], spacing=lc)
            net.add_boundary_pores(labels=['front', 'back'])

            Ps = net.pores('*boundary', mode='not')
            Ts = net.throats('*boundary', mode='not')
            BPs = net.pores('*boundary')
            BTs = net.throats('*boundary')
            geom = op.geometry.GenericGeometry(network=net, pores=Ps,
                                               throats=Ts)
            geom['pore.diameter'] = 2*rad
            geom['throat.diameter'] = 2*rad
            geom['pore.area'] = area
            geom['throat.area'] = area
            net['throat.centroid'] = throat_centroid(net)
            net['throat.vector'] = throat_vector(net)
            geom.add_model(propname='throat.ctc',
                           model=op.models.geometry.throat_length.ctc)
            geom['throat.length'] = geom['throat.ctc']*throat_portion
            geom.add_model(propname='throat.endpoints',
                           model=mods.straight_throat)
            geom.add_model(propname='throat.conduit_lengths',
                           model=gm.throat_length.conduit_lengths)
            boun = op.geometry.Boundary(network=net, pores=BPs, throats=BTs)
            air = op.phases.Air(network=net)
            air['pore.diffusivity'] = D_ab
            phys = op.physics.GenericPhysics(network=net, phase=air,
                                             geometry=geom)
            phys.add_model(propname='throat.diffusive_conductance',
                           model=pm.diffusive_conductance.ordinary_diffusion)
            physb = op.physics.GenericPhysics(network=net, phase=air,
                                              geometry=boun)
            physb.add_model(propname='throat.diffusive_conductance',
                            model=pm.diffusive_conductance.ordinary_diffusion)
            FD = op.algorithms.FickianDiffusion(network=net)
            FD.setup(phase=air)
            FD.set_value_BC(pores=net.pores('front_boundary'), values=1.0)
            FD.set_value_BC(pores=net.pores('back_boundary'), values=0.0)
            FD.run()
            D = FD.calc_effective_diffusivity(domain_area=area,
                                              domain_length=length)
            conds.append(D[0])
            nets.append(net)
        assert np.allclose(np.asarray(conds), D_ab)


if __name__ == '__main__':

    t = ThroatEndpointsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
