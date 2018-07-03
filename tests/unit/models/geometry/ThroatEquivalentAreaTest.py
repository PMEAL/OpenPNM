import openpnm as op
import openpnm.models.geometry.throat_equivalent_area as mods
import numpy as np
from numpy.testing import assert_approx_equal


class ThroatSurfaceAreaTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.air,
                                              geometry=self.geo)
        self.geo['pore.diameter'] = 0.2
        self.geo['throat.diameter'] = 0.1

    def test_spherical_pores(self):
        self.geo['throat.area'] = np.pi * self.geo['throat.diameter']**2 / 4
        self.geo.add_model(propname='throat.conduit_lengths',
                           model=op.models.geometry.throat_length.circular_pores)
        self.geo.add_model(propname='throat.equivalent_area',
                           model=mods.spherical_pores,
                           regen_mode='normal')
        AEqP1 = np.unique(self.geo['throat.equivalent_area.pore1'])
        AEqTr = np.unique(self.geo['throat.equivalent_area.throat'])
        AEqP2 = np.unique(self.geo['throat.equivalent_area.pore2'])
        assert_approx_equal(np.array([0.02065897]), AEqP1)
        assert_approx_equal(np.array([0.007853982]), AEqTr)
        assert_approx_equal(np.array([0.02065897]), AEqP2)

    def test_circular_pores(self):
        self.geo['throat.area'] = self.geo['throat.diameter']
        self.geo.add_model(propname='throat.conduit_lengths',
                           model=op.models.geometry.throat_length.circular_pores)
        self.geo.add_model(propname='throat.equivalent_area',
                           model=mods.circular_pores,
                           regen_mode='normal')
        AEqP1 = np.unique(self.geo['throat.equivalent_area.pore1'])
        AEqTr = np.unique(self.geo['throat.equivalent_area.throat'])
        AEqP2 = np.unique(self.geo['throat.equivalent_area.pore2'])
        assert_approx_equal(np.array([0.16539867]), AEqP1)
        assert_approx_equal(np.array([0.1]), AEqTr)
        assert_approx_equal(np.array([0.16539867]), AEqP2)


if __name__ == '__main__':

    t = ThroatSurfaceAreaTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
