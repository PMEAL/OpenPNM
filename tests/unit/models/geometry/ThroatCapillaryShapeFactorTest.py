import numpy as np
import openpnm as op


class ThroatCapillaryShapeFactorTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.air = op.phase.Air(network=self.net)
        self.net['throat.cross_sectional_area'] = 1.0
        self.net['throat.perimeter'] = np.pi
        self.net['throat.diameter'] = (4 / np.pi)**(0.5)

    def test_mason_morrow(self):
        mod = op.models.geometry.throat_capillary_shape_factor.mason_morrow
        self.net.add_model(propname='throat.capillary_shape_factor',
                           model=mod,
                           throat_perimeter='throat.perimeter',
                           throat_area='throat.cross_sectional_area',
                           regen_mode='normal')
        a = np.unique(self.net['throat.capillary_shape_factor'])
        b = np.array(0.10132118, ndmin=1)
        assert np.allclose(a, b)

    def test_jenkins_rao(self):
        mod = op.models.geometry.throat_capillary_shape_factor.jenkins_rao
        self.net.add_model(propname='throat.capillary_shape_factor',
                           model=mod,
                           throat_perimeter='throat.perimeter',
                           throat_area='throat.cross_sectional_area',
                           throat_diameter='throat.diameter',
                           regen_mode='normal')
        a = np.unique(self.net['throat.capillary_shape_factor'])
        b = np.array(0.88622693, ndmin=1)
        assert np.allclose(a, b)


if __name__ == '__main__':

    t = ThroatCapillaryShapeFactorTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
