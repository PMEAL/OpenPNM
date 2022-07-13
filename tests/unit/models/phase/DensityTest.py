import openpnm as op
from numpy.testing import assert_approx_equal


class DensityTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.phase = op.phase.Phase(network=self.net)
        self.phase['pore.temperature'] = 298.0  # K
        self.phase['pore.pressure'] = 101325.0  # Pa
        self.phase['pore.molecular_weight'] = 0.018  # kg/mol
        self.phase['pore.molar_density'] = 55539.0  # mol/m3
        self.phase['pore.salinity'] = 0.0  # ppt

    def test_standard(self):
        # Liquid water
        self.phase.add_model(propname='pore.density',
                             model=op.models.phase.density.standard)
        assert_approx_equal(self.phase['pore.density'].mean(), 999.702)

    def test_ideal_gas(self):
        # Water vapor
        self.phase.add_model(propname='pore.density',
                             model=op.models.phase.density.ideal_gas)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.density'].mean(), 0.73610248)

    def test_water(self):
        # Liquid water
        self.phase.add_model(propname='pore.density',
                             model=op.models.phase.density.water)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.density'].mean(), 996.9522)

    def teardown_class(self):
        del(self.phase)
        del(self.net)


if __name__ == '__main__':

    t = DensityTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
