import OpenPNM
import scipy as sp


class MolarDensityTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        self.phase = OpenPNM.Phases.GenericPhase(network=self.net)
        self.phase['pore.molecular_weight'] = 0.0291  # kg/mol
        self.phase['pore.density'] = 1.19  # kg/m3
        self.phase['pore.temperature'] = 298.0  # K
        self.phase['pore.pressure'] = 101325  # Pa
        self.phase['pore.critical_temperature'] = 132.65  # K
        self.phase['pore.critical_pressure'] = 3771000.0  # Pa

    def test_standard(self):
        f = OpenPNM.Phases.models.molar_density.standard
        self.phase.models.add(propname='pore.molar_density',
                              model=f)
        assert sp.allclose(self.phase['pore.molar_density'], 40.8934707)

    def test_ideal_gas(self):
        f = OpenPNM.Phases.models.molar_density.ideal_gas
        self.phase.models.add(propname='pore.molar_density',
                              model=f)
        assert sp.allclose(self.phase['pore.molar_density'], 40.8945824)

    def test_vanderwaals(self):
        f = OpenPNM.Phases.models.molar_density.vanderwaals
        self.phase.models.add(propname='pore.molar_density',
                              model=f)
        assert sp.allclose(self.phase['pore.molar_density'], 40.92524916)
