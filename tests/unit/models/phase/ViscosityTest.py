import openpnm as op
import scipy as sp
from numpy.testing import assert_approx_equal


class ViscosityTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.phase = op.phase.GenericPhase(network=self.net)
        self.phase['pore.temperature'] = 298.0  # K
        self.phase['pore.molecular_weight'] = 0.018  # kg/mol
        self.phase['pore.critical_temperature'] = 647.15  # K
        self.phase['pore.critical_volume'] = 1.805e-5  # m3/kmol
        self.phase['pore.salinity'] = 0  # g/kg

    def test_water(self):
        self.phase.add_model(propname='pore.viscosity',
                             model=op.models.phase.viscosity.water)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.viscosity'].mean(),
                            0.0008931909)

    def test_reynolds(self):
        self.phase.add_model(propname='pore.viscosity',
                             model=op.models.phase.viscosity.reynolds,
                             u0=0.001,
                             b=0.001)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.viscosity'].mean(),
                            0.001347161)

    def test_chung(self):
        self.phase.add_model(propname='pore.viscosity',
                             model=op.models.phase.viscosity.chung)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.viscosity'].mean(),
                            6.47289919e-05)


if __name__ == '__main__':

    t = ViscosityTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
