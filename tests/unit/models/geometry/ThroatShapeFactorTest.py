import numpy as np
import scipy as sp
import openpnm as op


class ThroatShapeFactorTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.air,
                                              geometry=self.geo)
        self.geo['throat.area'] = 1.0
        self.geo['throat.perimeter'] = np.pi
        self.geo['throat.diameter'] = (4/np.pi)**(1/2)

    def test_compactness(self):
        pass

    def test_mason_morrow(self):
        mod = op.models.geometry.throat_shape_factor.mason_morrow
        self.geo.add_model(propname='throat.shape_factor',
                           model=mod,
                           throat_perimeter='throat.perimeter',
                           throat_area='throat.area',
                           regen_mode='normal')
        a = np.unique(self.geo['throat.shape_factor'])
        b = np.array(0.10132118, ndmin=1)
        assert np.allclose(a, b)

    def test_jenkins_rao(self):
        mod = op.models.geometry.throat_shape_factor.jenkins_rao
        self.geo.add_model(propname='throat.shape_factor',
                           model=mod,
                           throat_perimeter='throat.perimeter',
                           throat_area='throat.area',
                           throat_diameter='throat.diameter',
                           regen_mode='normal')
        a = np.unique(self.geo['throat.shape_factor'])
        b = np.array(0.88622693, ndmin=1)
        assert np.allclose(a, b)


if __name__ == '__main__':

    t = ThroatShapeFactorTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
