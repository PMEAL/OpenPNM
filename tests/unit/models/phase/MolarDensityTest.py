import openpnm as op
from numpy.testing import assert_approx_equal


class MolarDensityTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.phase = op.phase.Phase(network=self.net)
        self.phase['pore.molecular_weight'] = 0.0291  # kg/mol
        self.phase['pore.density'] = 1.19  # kg/m3
        self.phase['pore.temperature'] = 298.0  # K
        self.phase['pore.pressure'] = 101325  # Pa
        self.phase['pore.critical_temperature'] = 132.65  # K
        self.phase['pore.critical_pressure'] = 3771000.0  # Pa

    def test_standard(self):
        f = op.models.phase.molar_density.standard
        self.phase.add_model(propname='pore.molar_density', model=f)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.molar_density'].mean(),
                            40.8934707)

    def test_ideal_gas(self):
        f = op.models.phase.molar_density.ideal_gas
        self.phase.add_model(propname='pore.molar_density', model=f)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.molar_density'].mean(),
                            40.8945824)

    def test_vanderwaals(self):
        f = op.models.phase.molar_density.vanderwaals
        self.phase.add_model(propname='pore.molar_density', model=f)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.molar_density'].mean(),
                            40.92524916)


if __name__ == '__main__':

    t = MolarDensityTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
