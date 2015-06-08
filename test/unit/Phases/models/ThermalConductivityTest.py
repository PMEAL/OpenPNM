import OpenPNM
import scipy as sp


class ThermalConductivityTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        self.phase = OpenPNM.Phases.GenericPhase(network=self.net)
        self.phase['pore.temperature'] = 298.0  # K
        self.phase['pore.salinity'] = 0.0  # g/kg
        self.phase['pore.viscosity'] = 0.001  # Pa.s
        self.phase['pore.critical_temperature'] = 647.15  # K
        self.phase['pore.molecular_weight'] = 0.018  # kg/mol
        self.phase['pore.boiling_point'] = 373.15  # K
        self.phase['pore.heat_capacity'] = 75.28  # J/mol K
        self.phase['pore.acentric_factor'] = 11.5  # J/mol K

    def test_water(self):
        f = OpenPNM.Phases.models.thermal_conductivity.water
        self.phase.models.add(propname='pore.thermal_conductivity',
                              model=f)
        assert sp.allclose(self.phase['pore.thermal_conductivity'], 0.61047611)

    def test_chung(self):
        f = OpenPNM.Phases.models.thermal_conductivity.chung
        self.phase.models.add(propname='pore.thermal_conductivity',
                              model=f)
        assert sp.allclose(self.phase['pore.thermal_conductivity'], 0.62063913)

    def test_sato(self):
        f = OpenPNM.Phases.models.thermal_conductivity.sato
        self.phase.models.add(propname='pore.thermal_conductivity',
                              model=f)
        assert sp.allclose(self.phase['pore.thermal_conductivity'], 0.29787023)
