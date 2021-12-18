import openpnm as op
import scipy as sp
from numpy.testing import assert_approx_equal


class ThermalConductivityTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.phase = op.phase.GenericPhase(network=self.net)
        self.phase['pore.temperature'] = 298.0  # K
        self.phase['pore.salinity'] = 0.0  # g/kg
        self.phase['pore.viscosity'] = 0.001  # Pa.s
        self.phase['pore.critical_temperature'] = 647.15  # K
        self.phase['pore.molecular_weight'] = 0.018  # kg/mol
        self.phase['pore.boiling_point'] = 373.15  # K
        self.phase['pore.heat_capacity'] = 75.28  # J/mol K
        self.phase['pore.acentric_factor'] = 11.5  # J/mol K

    def test_water(self):
        f = op.models.phase.thermal_conductivity.water
        self.phase.add_model(propname='pore.thermal_conductivity',
                             model=f)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.thermal_conductivity'].mean(),
                            0.61047611)

    def test_chung(self):
        f = op.models.phase.thermal_conductivity.chung
        self.phase.add_model(propname='pore.thermal_conductivity',
                             model=f)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.thermal_conductivity'].mean(),
                            0.62063913)

    def test_sato(self):
        f = op.models.phase.thermal_conductivity.sato
        self.phase.add_model(propname='pore.thermal_conductivity',
                             model=f)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.thermal_conductivity'].mean(),
                            0.29787023)


if __name__ == '__main__':

    t = ThermalConductivityTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
