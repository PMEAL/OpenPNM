import OpenPNM
import scipy as sp


class ViscosityTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        self.phase = OpenPNM.Phases.GenericPhase(network=self.net)
        self.phase['pore.temperature'] = 298.0  # K
        self.phase['pore.molecular_weight'] = 0.018  # kg/mol
        self.phase['pore.critical_temperature'] = 647.15  # K
        self.phase['pore.critical_volume'] = 1.805e-5  # m3/kmol
        self.phase['pore.salinity'] = 0  # g/kg

    def test_water(self):
        self.phase.models.add(propname='pore.viscosity',
                              model=OpenPNM.Phases.models.viscosity.water)
        assert sp.allclose(self.phase['pore.viscosity'], 0.00089319)

    def test_reynolds(self):
        self.phase.models.add(propname='pore.viscosity',
                              model=OpenPNM.Phases.models.viscosity.reynolds,
                              uo=0.001,
                              b=0.001)
        assert sp.allclose(self.phase['pore.viscosity'], 0.00134716)

    def test_chung(self):
        self.phase.models.add(propname='pore.viscosity',
                              model=OpenPNM.Phases.models.viscosity.chung)
        assert sp.allclose(self.phase['pore.viscosity'], 6.47289919e-05)
